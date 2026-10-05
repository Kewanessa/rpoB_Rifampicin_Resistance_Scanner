# Makefile for rpoB Rifampicin Resistance Scanner
# kewanessa
# CS50 Final Project
# AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
#
# Type 'make' to build, 'make clean' to remove compiled files.
# Type 'make test' to run the test suite.

CC = gcc
CFLAGS = -Wall -Wextra -Werror -g -O2 -std=c11
TARGET = amr_scan
ASAN_TARGET = amr_scan_asan

# Source files and their object files
SRCS = main.c fasta.c align.c translate.c variants.c mutdb.c report.c
HDRS = fasta.h align.h translate.h variants.h mutdb.h report.h
OBJS = $(SRCS:.c=.o)

# Default target: build the program
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

# Compile each .c to .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Dependencies (which headers each .c needs, including headers that
# those headers include: variants.h includes mutdb.h, report.h includes both)
main.o: main.c fasta.h align.h translate.h variants.h mutdb.h report.h
fasta.o: fasta.c fasta.h
align.o: align.c align.h
translate.o: translate.c translate.h
variants.o: variants.c variants.h mutdb.h translate.h
mutdb.o: mutdb.c mutdb.h
report.o: report.c report.h variants.h mutdb.h

# Clean up compiled files
clean:
	rm -f $(TARGET) $(ASAN_TARGET) $(OBJS)

# Run the test suite
test: $(TARGET)
	@echo "Running tests..."
	@bash tests/run_tests.sh

# Build a separate binary with AddressSanitizer, LeakSanitizer and
# UndefinedBehaviorSanitizer (the first error stops the program),
# then run the whole test suite on it
asan: $(ASAN_TARGET)
	@echo "Running tests with AddressSanitizer..."
	@PROG=./$(ASAN_TARGET) bash tests/run_tests.sh

$(ASAN_TARGET): $(SRCS) $(HDRS)
	$(CC) -Wall -Wextra -Werror -g -std=c11 -fsanitize=address,undefined \
		-fno-sanitize-recover=all -o $(ASAN_TARGET) $(SRCS) -lm

# Run with Valgrind to check for memory leaks
valgrind: $(TARGET)
	valgrind --leak-check=full --show-leak-kinds=all --error-exitcode=1 \
		./$(TARGET) tests/test_S450L.fasta "ncbi data/ncbi_dataset/data/gene.fna" mutations.tsv

.PHONY: clean test asan valgrind
