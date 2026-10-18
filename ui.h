/**
 * ui.h — Terminal output helpers
 *
 * Shared layout for the progress log and the report: numbered step
 * headings, aligned "label  value" fields, section headings with a
 * rule, word wrapping, and optional colour.
 *
 * Colour is used only when stdout is a terminal, and never when the
 * NO_COLOR environment variable is set (https://no-color.org) or
 * TERM=dumb, so redirected or piped output is always plain text.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: written with Anthropic Claude (Claude Code); see README.
 */

#ifndef UI_H
#define UI_H

// Lines are laid out to fit this many columns (an 80-column terminal with room to spare)
#define UI_WIDTH 72

// Let the compiler check printf-style format strings
#if defined(__GNUC__)
#define UI_PRINTF(fmt, args) __attribute__((format(printf, fmt, args)))
#else
#define UI_PRINTF(fmt, args)
#endif

// ANSI escape sequences, or "" when colour is off
typedef struct
{
    const char *reset;
    const char *bold;
    const char *dim;
    const char *red;
    const char *green;
    const char *yellow;
    const char *cyan;
} UiStyle;

extern UiStyle ui;

// Turn colour on if stdout is a terminal (call once, before printing).
void ui_init(void);

// Program name and version, underlined, at the top of the run.
void ui_banner(const char *title, const char *version);

// Heading for one pipeline step, e.g. "[1/6] Reading patient sequence".
void ui_step(int number, int total, const char *title);

// Detail line under a step heading: "File        tests/test_S450L.fasta".
void ui_detail(const char *label, const char *format, ...) UI_PRINTF(2, 3);

// Dimmed remark under a step heading, e.g. "this may take a moment".
void ui_note(const char *format, ...) UI_PRINTF(1, 2);

// Report section heading with a rule below it. A count >= 0 is shown at
// the right margin, in 'color' when non-zero.
void ui_heading(const char *title, int count, const char *color);

// Report field: "Sample:         patient_S450L".
void ui_field(const char *label, const char *format, ...) UI_PRINTF(2, 3);

// Horizontal line across the page, made of 'piece' (e.g. "─") in 'color'.
void ui_rule(const char *piece, const char *color);

/**
 * Print text from the current cursor column, breaking lines between words
 * so nothing passes UI_WIDTH, and wherever the text has a '\n'.
 * Continuation lines start at 'column'. Ends with a newline.
 */
void ui_wrap(const char *text, int column);

#endif // UI_H
