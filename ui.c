/**
 * ui.c — Terminal output helpers
 *
 * THE LAYOUT GRID (columns):
 *
 *   2   step headings, section headings, rules, the verdict
 *   4   report fields and variant bullets
 *   6   lines under a variant (interpretation, source)
 *   8   step details, lined up under the step title
 *
 * Lines are laid out for UI_WIDTH columns; long values wrap under themselves.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: written with Anthropic Claude (Claude Code); see README.
 */

// isatty() and fileno() are POSIX, not part of strict C11
#define _POSIX_C_SOURCE 200809L

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ui.h"

// Step details: indent (under the "[1/6] " heading) and label width
#define STEP_INDENT       8
#define STEP_LABEL_WIDTH  12

// Report fields: indent and label width ("Sequenced Span: " is the longest)
#define FIELD_INDENT      4
#define FIELD_LABEL_WIDTH 16

// No colour until ui_init() decides the terminal can show it
UiStyle ui = {"", "", "", "", "", "", ""};

void ui_init(void)
{
    const char *no_color = getenv("NO_COLOR");
    const char *term = getenv("TERM");

    if (!isatty(fileno(stdout)) ||
        (no_color != NULL && no_color[0] != '\0') ||
        (term != NULL && strcmp(term, "dumb") == 0))
    {
        return;
    }

    ui.reset  = "\033[0m";
    ui.bold   = "\033[1m";
    ui.dim    = "\033[2m";
    ui.red    = "\033[31m";
    ui.green  = "\033[32m";
    ui.yellow = "\033[33m";
    ui.cyan   = "\033[36m";
}

/**
 * Number of screen columns taken by the first len bytes of s.
 * ANSI colour codes take none, and a multi-byte UTF-8 character
 * (such as ✓) takes one: only its first byte is counted.
 */
static int visible_width(const char *s, int len)
{
    int width = 0;
    for (int i = 0; i < len; i++)
    {
        if (s[i] == '\033')
        {
            while (i < len && s[i] != 'm') i++;
        }
        else if (((unsigned char) s[i] & 0xC0) != 0x80)
        {
            width++;
        }
    }
    return width;
}

void ui_wrap(const char *text, int column)
{
    int cursor = column;  // Screen column the next character lands in
    int line_empty = 1;   // No word printed on this line yet
    const char *p = text;

    while (*p != '\0')
    {
        // Next word: everything up to the next space or newline
        while (*p == ' ') p++;
        if (*p == '\0') break;
        int len = strcspn(p, " \n");
        int width = visible_width(p, len);

        // Move to a new line if the text asks for one, or if the word would
        // pass the right margin (a word longer than a whole line is printed anyway)
        if (*p == '\n')
        {
            printf("\n%*s", column, "");
            cursor = column;
            line_empty = 1;
            p++;
            continue;
        }
        if (!line_empty && cursor + 1 + width > UI_WIDTH)
        {
            printf("\n%*s", column, "");
            cursor = column;
            line_empty = 1;
        }
        if (!line_empty)
        {
            printf(" ");
            cursor++;
        }
        printf("%.*s", len, p);
        cursor += width;
        line_empty = 0;
        p += len;
    }
    printf("\n");
}

/**
 * Print a dimmed label padded to label_width, then the formatted value,
 * wrapped so that continuation lines line up with the value.
 */
static void print_field(int indent, int label_width, const char *label,
                        const char *format, va_list args) UI_PRINTF(4, 0);

static void print_field(int indent, int label_width, const char *label,
                        const char *format, va_list args)
{
    printf("%*s%s%-*s%s", indent, "", ui.dim, label_width, label, ui.reset);

    // A label wider than label_width pushes the value right
    int column = indent + label_width;
    int label_len = strlen(label);
    if (label_len > label_width)
    {
        column = indent + label_len;
    }

    // Format the value into a buffer of the right size, then wrap it
    va_list copy;
    va_copy(copy, args);
    int len = vsnprintf(NULL, 0, format, copy);
    va_end(copy);

    char *text = (len >= 0) ? malloc(len + 1) : NULL;
    if (text == NULL)
    {
        // Out of memory: print the value as it is, without wrapping
        vprintf(format, args);
        printf("\n");
        return;
    }
    vsnprintf(text, len + 1, format, args);
    ui_wrap(text, column);
    free(text);
}

void ui_banner(const char *title, const char *version)
{
    printf("\n  %s%s%s  %sv%s%s\n", ui.bold, title, ui.reset, ui.dim, version, ui.reset);
    ui_rule("─", ui.dim);
}

void ui_step(int number, int total, const char *title)
{
    printf("\n  %s%s[%d/%d]%s %s%s%s\n", ui.bold, ui.cyan, number, total, ui.reset,
           ui.bold, title, ui.reset);

    // Show progress right away, even when stdout is a pipe or a file
    fflush(stdout);
}

void ui_detail(const char *label, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    print_field(STEP_INDENT, STEP_LABEL_WIDTH, label, format, args);
    va_end(args);
    fflush(stdout);
}

void ui_note(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    printf("%*s%s", STEP_INDENT, "", ui.dim);
    vprintf(format, args);
    printf("%s\n", ui.reset);
    va_end(args);
    fflush(stdout);
}

void ui_heading(const char *title, int count, const char *color)
{
    printf("\n  %s%s%s", ui.bold, title, ui.reset);
    if (count >= 0)
    {
        // Right-align the count with the end of the rule below
        int space = UI_WIDTH - 2 - (int) strlen(title);
        if (count > 0)
        {
            printf("%s%s%*d%s", ui.bold, color, space, count, ui.reset);
        }
        else
        {
            printf("%s%*d%s", ui.dim, space, count, ui.reset);
        }
    }
    printf("\n");
    ui_rule("─", ui.dim);
}

void ui_field(const char *label, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    print_field(FIELD_INDENT, FIELD_LABEL_WIDTH, label, format, args);
    va_end(args);
}

void ui_rule(const char *piece, const char *color)
{
    printf("  %s", color);
    for (int i = 2; i < UI_WIDTH; i++)
    {
        printf("%s", piece);
    }
    printf("%s\n", ui.reset);
}
