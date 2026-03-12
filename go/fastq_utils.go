package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// FastqRecord represents a FASTQ record
type FastqRecord struct {
	Name     string
	Sequence string
	Plus     string
	Quality  string
}

// FastqReader reads FASTQ files
type FastqReader struct {
	reader *bufio.Reader
}

// NewFastqReader creates a new FASTQ reader
func NewFastqReader(filename string) (*FastqReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	var r io.Reader = f
	if strings.HasSuffix(filename, ".gz") {
		gr, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		r = gr
	}

	return &FastqReader{reader: bufio.NewReader(r)}, nil
}

// Read reads the next FASTQ record
func (r *FastqReader) Read() (*FastqRecord, error) {
	// Read name line (starts with @)
	nameLine, err := r.reader.ReadString('\n')
	if err != nil {
		return nil, err
	}
	nameLine = strings.TrimSpace(nameLine)
	if !strings.HasPrefix(nameLine, "@") {
		return nil, fmt.Errorf("invalid FASTQ: expected @ at start of record")
	}
	name := nameLine[1:]

	// Read sequence
	sequence, err := r.reader.ReadString('\n')
	if err != nil {
		return nil, err
	}
	sequence = strings.TrimSpace(sequence)

	// Read plus line
	plus, err := r.reader.ReadString('\n')
	if err != nil {
		return nil, err
	}
	plus = strings.TrimSpace(plus)

	// Read quality
	quality, err := r.reader.ReadString('\n')
	if err != nil {
		return nil, err
	}
	quality = strings.TrimSpace(quality)

	return &FastqRecord{
		Name:     name,
		Sequence: sequence,
		Plus:     plus,
		Quality:  quality,
	}, nil
}

// FastqWriter writes FASTQ records
type FastqWriter struct {
	writer *bufio.Writer
}

// NewFastqWriter creates a new FASTQ writer
func NewFastqWriter(filename string, compress bool) (*FastqWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}

	var w io.Writer = f
	if compress {
		gw := gzip.NewWriter(f)
		w = gw
	}

	return &FastqWriter{writer: bufio.NewWriter(w)}, nil
}

// Write writes a FASTQ record
func (w *FastqWriter) Write(rec *FastqRecord) error {
	_, err := fmt.Fprintf(w.writer, "@%s\n%s\n%s\n%s\n",
		rec.Name, rec.Sequence, rec.Plus, rec.Quality)
	return err
}

// Flush flushes the writer
func (w *FastqWriter) Flush() error {
	return w.writer.Flush()
}

// FilterFastqByLength filters FASTQ records by length
func FilterFastqByLength(inputPath, outputPath string, minLen, maxLen int, compress bool) error {
	reader, err := NewFastqReader(inputPath)
	if err != nil {
		return err
	}

	writer, err := NewFastqWriter(outputPath, compress)
	if err != nil {
		return err
	}
	defer writer.Flush()

	count := 0
	written := 0
	for {
		rec, err := reader.Read()
		if err != nil {
			break
		}
		count++

		l := len(rec.Sequence)
		if l >= minLen && (maxLen == 0 || l <= maxLen) {
			writer.Write(rec)
			written++
		}
	}

	fmt.Printf("Processed %d records, wrote %d\n", count, written)
	return nil
}
