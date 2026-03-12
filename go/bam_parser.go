package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// ReadRecord represents a single read with metadata
type ReadRecord struct {
	ReadID          string
	Sequence        string
	Quality         string
	Barcode         string
	UMI             string
	MapQ            int
	IsUnmapped      bool
	IsSecondary     bool
	IsSupplementary bool
	ReferenceName   string
	ReferenceStart  int
}

// ToFastq converts record to FASTQ format
func (r *ReadRecord) ToFastq() string {
	return fmt.Sprintf("@%s\n%s\n+\n%s\n", r.ReadID, r.Sequence, r.Quality)
}

// ToFasta converts record to FASTA format
func (r *ReadRecord) ToFasta() string {
	return fmt.Sprintf(">%s\n%s\n", r.ReadID, r.Sequence)
}

// BAMParser parses scRNA-seq BAM files
type BAMParser struct {
	bamPath     string
	barcodeTag  string
	umiTag      string
	chunkSize   int
	numWorkers  int
}

// NewBAMParser creates a new BAM parser
func NewBAMParser(bamPath, barcodeTag, umiTag string, chunkSize, numWorkers int) *BAMParser {
	if numWorkers <= 0 {
		numWorkers = runtime.NumCPU()
	}
	return &BAMParser{
		bamPath:    bamPath,
		barcodeTag: barcodeTag,
		umiTag:     umiTag,
		chunkSize:  chunkSize,
		numWorkers: numWorkers,
	}
}

// extractBarcode extracts cell barcode from SAM record
func (p *BAMParser) extractBarcode(rec *sam.Record) string {
	for _, aux := range rec.AuxFields {
		tag := aux.Tag()
		if string(tag[:]) == p.barcodeTag {
			val := strings.TrimPrefix(aux.String(), p.barcodeTag+":Z:")
			// Remove -1 suffix if present
			if idx := strings.Index(val, "-"); idx != -1 {
				val = val[:idx]
			}
			return val
		}
	}
	return ""
}

// extractUMI extracts UMI from SAM record
func (p *BAMParser) extractUMI(rec *sam.Record) string {
	for _, aux := range rec.AuxFields {
		tag := aux.Tag()
		if string(tag[:]) == p.umiTag {
			return strings.TrimPrefix(aux.String(), p.umiTag+":Z:")
		}
	}
	return ""
}

// recordToStruct converts SAM record to ReadRecord
func (p *BAMParser) recordToStruct(rec *sam.Record) *ReadRecord {
	refName := ""
	refStart := 0
	if rec.Ref != nil {
		refName = rec.Ref.Name()
		refStart = rec.Start()
	}

	seq := ""
	if rec.Seq.Length > 0 {
		seqBytes := rec.Seq.Expand()
		seq = string(seqBytes)
	}

	qual := ""
	if len(rec.Qual) > 0 {
		qual = string(rec.Qual)
	}

	return &ReadRecord{
		ReadID:          rec.Name,
		Sequence:        seq,
		Quality:         qual,
		Barcode:         p.extractBarcode(rec),
		UMI:             p.extractUMI(rec),
		MapQ:            int(rec.MapQ),
		IsUnmapped:      rec.Flags&sam.Unmapped != 0,
		IsSecondary:     rec.Flags&sam.Secondary != 0,
		IsSupplementary: rec.Flags&sam.Supplementary != 0,
		ReferenceName:   refName,
		ReferenceStart:  refStart,
	}
}

// Stats holds parsing statistics
type Stats struct {
	TotalReads       int64
	MappedReads      int64
	UnmappedReads    int64
	LowMapQReads     int64
	ReadsWithBarcode int64
	ReadsWithUMI     int64
	CandidateReads   int64
}

// ExtractCandidateReads extracts candidate microbial reads
func (p *BAMParser) ExtractCandidateReads(mapqThreshold int, includeUnmapped, includeLowMapQ bool) ([]*ReadRecord, *Stats, error) {
	f, err := os.Open(p.bamPath)
	if err != nil {
		return nil, nil, fmt.Errorf("opening BAM: %w", err)
	}
	defer f.Close()

	br, err := bam.NewReader(f, p.numWorkers)
	if err != nil {
		return nil, nil, fmt.Errorf("creating BAM reader: %w", err)
	}
	defer br.Close()

	stats := &Stats{}
	var candidates []*ReadRecord

	for {
		rec, err := br.Read()
		if err != nil {
			break
		}

		stats.TotalReads++

		// Skip secondary and supplementary alignments
		if rec.Flags&sam.Secondary != 0 || rec.Flags&sam.Supplementary != 0 {
			continue
		}

		// Track statistics
		if rec.Flags&sam.Unmapped != 0 {
			stats.UnmappedReads++
		} else {
			stats.MappedReads++
			if int(rec.MapQ) < mapqThreshold {
				stats.LowMapQReads++
			}
		}

		// Check barcode and UMI
		if p.extractBarcode(rec) != "" {
			stats.ReadsWithBarcode++
		}
		if p.extractUMI(rec) != "" {
			stats.ReadsWithUMI++
		}

		// Determine if candidate
		isCandidate := false
		if includeUnmapped && rec.Flags&sam.Unmapped != 0 {
			isCandidate = true
		} else if includeLowMapQ && rec.Flags&sam.Unmapped == 0 && int(rec.MapQ) < mapqThreshold {
			isCandidate = true
		}

		if isCandidate {
			stats.CandidateReads++
			candidates = append(candidates, p.recordToStruct(rec))
		}
	}

	return candidates, stats, nil
}

// ExtractCandidateReadsParallel extracts candidates using multiple workers
func (p *BAMParser) ExtractCandidateReadsParallel(mapqThreshold int, includeUnmapped, includeLowMapQ bool) ([]*ReadRecord, *Stats, error) {
	f, err := os.Open(p.bamPath)
	if err != nil {
		return nil, nil, fmt.Errorf("opening BAM: %w", err)
	}
	defer f.Close()

	br, err := bam.NewReader(f, 0)
	if err != nil {
		return nil, nil, fmt.Errorf("creating BAM reader: %w", err)
	}
	defer br.Close()

	// Channels for parallel processing
	recordChan := make(chan *sam.Record, p.chunkSize*2)
	resultChan := make(chan *ReadRecord, p.chunkSize*2)
	var wg sync.WaitGroup

	// Start workers
	for i := 0; i < p.numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for rec := range recordChan {
				// Skip secondary and supplementary
				if rec.Flags&sam.Secondary != 0 || rec.Flags&sam.Supplementary != 0 {
					continue
				}

				// Check if candidate
				isCandidate := false
				if includeUnmapped && rec.Flags&sam.Unmapped != 0 {
					isCandidate = true
				} else if includeLowMapQ && rec.Flags&sam.Unmapped == 0 && int(rec.MapQ) < mapqThreshold {
					isCandidate = true
				}

				if isCandidate {
					resultChan <- p.recordToStruct(rec)
				}
			}
		}()
	}

	// Collector goroutine
	var candidates []*ReadRecord
	var collectWg sync.WaitGroup
	collectWg.Add(1)
	go func() {
		defer collectWg.Done()
		for rec := range resultChan {
			candidates = append(candidates, rec)
		}
	}()

	// Read records and send to workers
	stats := &Stats{}
	for {
		rec, err := br.Read()
		if err != nil {
			break
		}

		stats.TotalReads++

		if rec.Flags&sam.Unmapped != 0 {
			stats.UnmappedReads++
		} else {
			stats.MappedReads++
			if int(rec.MapQ) < mapqThreshold {
				stats.LowMapQReads++
			}
		}

		if p.extractBarcode(rec) != "" {
			stats.ReadsWithBarcode++
		}
		if p.extractUMI(rec) != "" {
			stats.ReadsWithUMI++
		}

		recordChan <- rec
	}

	close(recordChan)
	wg.Wait()
	close(resultChan)
	collectWg.Wait()

	stats.CandidateReads = int64(len(candidates))

	return candidates, stats, nil
}

// WriteToFastq writes records to FASTQ file
func WriteToFastq(records []*ReadRecord, outputPath string, compress bool) error {
	f, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer f.Close()

	var w *bufio.Writer
	if compress {
		gw := gzip.NewWriter(f)
		defer gw.Close()
		w = bufio.NewWriter(gw)
	} else {
		w = bufio.NewWriter(f)
	}
	defer w.Flush()

	for _, rec := range records {
		w.WriteString(rec.ToFastq())
	}

	return nil
}

// WriteToFasta writes records to FASTA file
func WriteToFasta(records []*ReadRecord, outputPath string, compress bool) error {
	f, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer f.Close()

	var w *bufio.Writer
	if compress {
		gw := gzip.NewWriter(f)
		defer gw.Close()
		w = bufio.NewWriter(gw)
	} else {
		w = bufio.NewWriter(f)
	}
	defer w.Flush()

	for _, rec := range records {
		w.WriteString(rec.ToFasta())
	}

	return nil
}

// WriteBarcodes writes unique barcodes to file
func WriteBarcodes(records []*ReadRecord, outputPath string) error {
	barcodes := make(map[string]bool)
	for _, rec := range records {
		if rec.Barcode != "" {
			barcodes[rec.Barcode] = true
		}
	}

	f, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()

	for barcode := range barcodes {
		w.Write([]string{barcode})
	}

	return nil
}

func main() {
	var (
		bamPath        = flag.String("i", "", "Input BAM file")
		outputFastq    = flag.String("o-fastq", "", "Output FASTQ file")
		outputFasta    = flag.String("o-fasta", "", "Output FASTA file")
		outputBarcodes = flag.String("o-barcodes", "", "Output barcodes file")
		barcodeTag     = flag.String("barcode-tag", "CB", "Cell barcode tag")
		umiTag         = flag.String("umi-tag", "UB", "UMI tag")
		mapqThreshold  = flag.Int("mapq", 30, "MAPQ threshold")
		chunkSize      = flag.Int("chunk-size", 10000, "Chunk size")
		numWorkers     = flag.Int("threads", 0, "Number of workers (0 = auto)")
		parallel       = flag.Bool("parallel", false, "Use parallel processing")
		compress       = flag.Bool("compress", false, "Compress output")
	)
	flag.Parse()

	if *bamPath == "" {
		log.Fatal("Please provide input BAM file with -i")
	}

	parser := NewBAMParser(*bamPath, *barcodeTag, *umiTag, *chunkSize, *numWorkers)

	log.Printf("Parsing BAM file: %s", *bamPath)

	var candidates []*ReadRecord
	var stats *Stats
	var err error

	if *parallel {
		candidates, stats, err = parser.ExtractCandidateReadsParallel(*mapqThreshold, true, true)
	} else {
		candidates, stats, err = parser.ExtractCandidateReads(*mapqThreshold, true, true)
	}

	if err != nil {
		log.Fatalf("Error parsing BAM: %v", err)
	}

	log.Printf("Parsing complete:")
	log.Printf("  Total reads: %d", stats.TotalReads)
	log.Printf("  Mapped reads: %d", stats.MappedReads)
	log.Printf("  Unmapped reads: %d", stats.UnmappedReads)
	log.Printf("  Low MAPQ reads: %d", stats.LowMapQReads)
	log.Printf("  Reads with barcode: %d", stats.ReadsWithBarcode)
	log.Printf("  Reads with UMI: %d", stats.ReadsWithUMI)
	log.Printf("  Candidate reads: %d", stats.CandidateReads)

	// Write outputs
	if *outputFastq != "" {
		log.Printf("Writing FASTQ: %s", *outputFastq)
		if err := WriteToFastq(candidates, *outputFastq, *compress); err != nil {
			log.Fatalf("Error writing FASTQ: %v", err)
		}
	}

	if *outputFasta != "" {
		log.Printf("Writing FASTA: %s", *outputFasta)
		if err := WriteToFasta(candidates, *outputFasta, *compress); err != nil {
			log.Fatalf("Error writing FASTA: %v", err)
		}
	}

	if *outputBarcodes != "" {
		log.Printf("Writing barcodes: %s", *outputBarcodes)
		if err := WriteBarcodes(candidates, *outputBarcodes); err != nil {
			log.Fatalf("Error writing barcodes: %v", err)
		}
	}

	log.Println("Done!")
}
