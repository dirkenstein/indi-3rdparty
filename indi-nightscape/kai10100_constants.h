#ifndef KAI10100_CONSTANTS_H
#define KAI10100_CONSTANTS_H

#define KAI10100_MAX_X 3868
#define KAI10100_MAX_Y 2892
#define KAI10100_ACTIVE_X 3760

// 4 dummy + 8 active buf
#define KAI10100_X_PREAMBLE (4 + 8)
// 8 active buf + 2 dark dummy + 74 dark + 12 dark dummy
#define KAI10100_X_POSTAMBLE (8 + 2 + 74 + 12)

// 4 gray + 8 active buffer
#define KAI10100_Y_PREAMBLE (4 + 8)
// 8 active buffer + 4 gray + 2 dark dummy + 20 dark + 2 dark dummy + 4 gray
#define KAI10100_Y_POSTAMBLE (8 + 4 + 2 + 20 + 2 + 4)

#define KAI10100_ACTIVE_Y 2840
#define KAI10100_HALF_Y   1428
#define KAI10100_QUARTER_Y 718
#define KAI10100_INTERLACE 716
#endif