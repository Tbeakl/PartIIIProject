/*
	This file is part of the ChipWhisperer Example Targets
	Copyright (C) 2012-2017 NewAE Technology Inc.

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hal.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "simpleserial.h"

// ChaCha implementation taken from https://github.com/rweather/lwc-finalists/blob/5d2b22c9ff7744be429cabda0c078ea5b7b6f79e/src/individual/ASCON_masked/aead-random.c
#define RO12(state, loc) \
	do                   \
	{                    \
		RO8(state, loc); \
		RO4(state, loc); \
	} while (0)

#define RO4(state, loc)                                                     \
	do                                                                      \
	{                                                                       \
		uint8_t _temp = state[4 * loc];                                     \
		state[4 * loc] = (state[4 * loc] << 4) | (state[4 * loc + 3] >> 4); \
		uint8_t _temp2 = state[4 * loc + 1];                                \
		state[4 * loc + 1] = (state[4 * loc + 1] << 4) | (_temp >> 4);      \
		_temp = state[4 * loc + 2];                                         \
		state[4 * loc + 2] = (state[4 * loc + 2] << 4) | (_temp2 >> 4);     \
		state[4 * loc + 3] = (state[4 * loc + 3] << 4) | (_temp >> 4);      \
	} while (0)

#define RO7(state, loc)                                                     \
	do                                                                      \
	{                                                                       \
		uint8_t _temp = state[4 * loc];                                     \
		state[4 * loc] = (state[4 * loc] << 7) | (state[4 * loc + 3] >> 1); \
		uint8_t _temp2 = state[4 * loc + 1];                                \
		state[4 * loc + 1] = (state[4 * loc + 1] << 7) | (_temp >> 1);      \
		_temp = state[4 * loc + 2];                                         \
		state[4 * loc + 2] = (state[4 * loc + 2] << 7) | (_temp2 >> 1);     \
		state[4 * loc + 3] = (state[4 * loc + 3] << 7) | (_temp >> 1);      \
	} while (0)

#define RO8(state, loc)                      \
	do                                       \
	{                                        \
		uint8_t _temp = state[4 * loc + 1];  \
		state[4 * loc + 1] = state[4 * loc]; \
		uint8_t _temp2 = state[4 * loc + 2]; \
		state[4 * loc + 2] = _temp;          \
		_temp = state[4 * loc + 3];          \
		state[4 * loc + 3] = _temp2;         \
		state[4 * loc] = _temp;              \
	} while (0)

#define RO16(state, loc)                         \
	do                                           \
	{                                            \
		uint8_t _temp = state[4 * loc];          \
		state[4 * loc] = state[4 * loc + 2];     \
		state[4 * loc + 2] = _temp;              \
		_temp = state[4 * loc + 1];              \
		state[4 * loc + 1] = state[4 * loc + 3]; \
		state[4 * loc + 3] = _temp;              \
	} while (0)

// #define ADD(state, output_loc, input_b_loc)                      \
// 	do                                                           \
// 	{                                                            \
// 		__asm("ADD %[out_a_1],%[in_b_1],%[in_a_1]\n\t"           \
// 			  "ADDS %[out_a_1], %[out_a_1], #4294967040\n\t"     \
// 			  "ADC %[out_a_2],%[in_b_2],%[in_a_2]\n\t"           \
// 			  "ADDS %[out_a_2], %[out_a_2], #4294967040\n\t"     \
// 			  "ADC %[out_a_3],%[in_b_3],%[in_a_3]\n\t"           \
// 			  "ADDS %[out_a_3], %[out_a_3], #4294967040\n\t"     \
// 			  "ADC %[out_a_4],%[in_b_4],%[in_a_4]\n\t"           \
// 			  "AND %[out_a_1], %[out_a_1], #255\n\t"             \
// 			  "AND %[out_a_2], %[out_a_2], #255\n\t"             \
// 			  "AND %[out_a_3], %[out_a_3], #255\n\t"             \
// 			  "AND %[out_a_4], %[out_a_4], #255\n\t"             \
// 			  : [out_a_1] "=r"(state[4 * output_loc]),           \
// 				[out_a_2] "=r"(state[4 * output_loc + 1]),       \
// 				[out_a_3] "=r"(state[4 * output_loc + 2]),       \
// 				[out_a_4] "=r"(state[4 * output_loc + 3])        \
// 			  :                                                  \
// 			  [in_a_1] "r"(state[4 * output_loc]),               \
// 			  [in_a_2] "r"(state[4 * output_loc + 1]),           \
// 			  [in_a_3] "r"(state[4 * output_loc + 2]),           \
// 			  [in_a_4] "r"(state[4 * output_loc + 3]),           \
// 			  [in_b_1] "r"(state[4 * input_b_loc]),      \
// 			  [in_b_2] "r"(state[4 * input_b_loc + 1]),  \
// 			  [in_b_3] "r"(state[4 * input_b_loc + 2]),  \
// 			  [in_b_4] "r"(state[4 * input_b_loc + 3])); \
// 	} while (0)

// #define ADD(state, initial_state, output_loc, input_b_loc)       \
// 	do                                                           \
// 	{                                                            \
// 		__asm("ADD %[out_a_1],%[in_b_1],%[in_a_1]\n\t"           \
// 			  "LSLS %[out_a_1], %[out_a_1], #24\n\t"             \
// 			  "ADC %[out_a_2],%[in_b_2],%[in_a_2]\n\t"           \
// 			  "LSLS %[out_a_2], %[out_a_2], #24\n\t"             \
// 			  "ADC %[out_a_3],%[in_b_3],%[in_a_3]\n\t"           \
// 			  "LSLS %[out_a_3], %[out_a_3], #24\n\t"             \
// 			  "ADC %[out_a_4],%[in_b_4],%[in_a_4]\n\t"           \
// 			  "LSR %[out_a_1], %[out_a_1], #24\n\t"              \
// 			  "LSR %[out_a_2], %[out_a_2], #24\n\t"              \
// 			  "LSR %[out_a_3], %[out_a_3], #24\n\t"              \
// 			  "AND %[out_a_4], %[out_a_4], #255\n\t"             \
// 			  : [out_a_1] "=r"(state[4 * output_loc]),           \
// 				[out_a_2] "=r"(state[4 * output_loc + 1]),       \
// 				[out_a_3] "=r"(state[4 * output_loc + 2]),       \
// 				[out_a_4] "=r"(state[4 * output_loc + 3])        \
// 			  :                                                  \
// 			  [in_a_1] "r"(initial_state[4 * output_loc]),       \
// 			  [in_a_2] "r"(initial_state[4 * output_loc + 1]),   \
// 			  [in_a_3] "r"(initial_state[4 * output_loc + 2]),   \
// 			  [in_a_4] "r"(initial_state[4 * output_loc + 3]),   \
// 			  [in_b_1] "r"(initial_state[4 * input_b_loc]),      \
// 			  [in_b_2] "r"(initial_state[4 * input_b_loc + 1]),  \
// 			  [in_b_3] "r"(initial_state[4 * input_b_loc + 2]),  \
// 			  [in_b_4] "r"(initial_state[4 * input_b_loc + 3])); \
// 	} while (0)

// #define ADD_END(state, initial_state, output_loc, input_b_loc)   \
// 	do                                                           \
// 	{                                                            \
// 		__asm("ADD %[out_a_1],%[in_b_1],%[in_a_1]\n\t"           \
// 			  "LSLS %[out_a_1], %[out_a_1], #24\n\t"             \
// 			  "ADC %[out_a_2],%[in_b_2],%[in_a_2]\n\t"           \
// 			  "LSLS %[out_a_2], %[out_a_2], #24\n\t"             \
// 			  "ADC %[out_a_3],%[in_b_3],%[in_a_3]\n\t"           \
// 			  "LSLS %[out_a_3], %[out_a_3], #24\n\t"             \
// 			  "ADC %[out_a_4],%[in_b_4],%[in_a_4]\n\t"           \
// 			  "LSR %[out_a_1], %[out_a_1], #24\n\t"              \
// 			  "LSR %[out_a_2], %[out_a_2], #24\n\t"              \
// 			  "LSR %[out_a_3], %[out_a_3], #24\n\t"              \
// 			  "AND %[out_a_4], %[out_a_4], #255\n\t"             \
// 			  : [out_a_1] "=r"(state[4 * output_loc]),           \
// 				[out_a_2] "=r"(state[4 * output_loc + 1]),       \
// 				[out_a_3] "=r"(state[4 * output_loc + 2]),       \
// 				[out_a_4] "=r"(state[4 * output_loc + 3])        \
// 			  :                                                  \
// 			  [in_a_1] "r"(state[4 * output_loc]),               \
// 			  [in_a_2] "r"(state[4 * output_loc + 1]),           \
// 			  [in_a_3] "r"(state[4 * output_loc + 2]),           \
// 			  [in_a_4] "r"(state[4 * output_loc + 3]),           \
// 			  [in_b_1] "r"(initial_state[4 * input_b_loc]),      \
// 			  [in_b_2] "r"(initial_state[4 * input_b_loc + 1]),  \
// 			  [in_b_3] "r"(initial_state[4 * input_b_loc + 2]),  \
// 			  [in_b_4] "r"(initial_state[4 * input_b_loc + 3])); \
// 	} while (0)

#define ADD(state, initial_state, output_loc, input_b_loc)      \
	do                                                          \
	{                                                           \
		uint32_t a_1 = initial_state[4 * output_loc];           \
		uint32_t a_2 = initial_state[4 * output_loc + 1];       \
		uint32_t a_3 = initial_state[4 * output_loc + 2];       \
		uint32_t a_4 = initial_state[4 * output_loc + 3];       \
		a_1 += initial_state[4 * input_b_loc];                  \
		uint8_t carry_flag = a_1 >> 8;                          \
		a_2 += carry_flag + initial_state[4 * input_b_loc + 1]; \
		carry_flag = a_2 >> 8;                                  \
		a_3 += carry_flag + initial_state[4 * input_b_loc + 2]; \
		carry_flag = a_3 >> 8;                                  \
		a_4 += carry_flag + initial_state[4 * input_b_loc + 3]; \
		state[4 * output_loc] = a_1 & 0xFF;                     \
		state[4 * output_loc + 1] = a_2 & 0xFF;                 \
		state[4 * output_loc + 2] = a_3 & 0xFF;                 \
		state[4 * output_loc + 3] = a_4 & 0xFF;                 \
	} while (0)

#define ADD_END(state, initial_state, output_loc, input_b_loc)  \
	do                                                          \
	{                                                           \
		uint32_t a_1 = state[4 * output_loc];                   \
		uint32_t a_2 = state[4 * output_loc + 1];               \
		uint32_t a_3 = state[4 * output_loc + 2];               \
		uint32_t a_4 = state[4 * output_loc + 3];               \
		a_1 += initial_state[4 * input_b_loc];                  \
		uint8_t carry_flag = a_1 >> 8;                          \
		a_2 += carry_flag + initial_state[4 * input_b_loc + 1]; \
		carry_flag = a_2 >> 8;                                  \
		a_3 += carry_flag + initial_state[4 * input_b_loc + 2]; \
		carry_flag = a_3 >> 8;                                  \
		a_4 += carry_flag + initial_state[4 * input_b_loc + 3]; \
		state[4 * output_loc] = a_1 & 0xFF;                     \
		state[4 * output_loc + 1] = a_2 & 0xFF;                 \
		state[4 * output_loc + 2] = a_3 & 0xFF;                 \
		state[4 * output_loc + 3] = a_4 & 0xFF;                 \
	} while (0)

#define XOR(state, output_loc, input_a_loc, input_b_loc)                                     \
	do                                                                                       \
	{                                                                                        \
		state[4 * output_loc] = state[4 * input_a_loc] ^ state[4 * input_b_loc];             \
		state[4 * output_loc + 1] = state[4 * input_a_loc + 1] ^ state[4 * input_b_loc + 1]; \
		state[4 * output_loc + 2] = state[4 * input_a_loc + 2] ^ state[4 * input_b_loc + 2]; \
		state[4 * output_loc + 3] = state[4 * input_a_loc + 3] ^ state[4 * input_b_loc + 3]; \
	} while (0)

/* Perform a ChaCha quarter round operation */
#define quarterRound_indices(state, a, b, c, d) \
	do                                          \
	{                                           \
		ADD(state, state, a, b);                \
		XOR(state, d, d, a);                    \
		RO16(state, d);                         \
		ADD(state, state, c, d);                \
		XOR(state, b, b, c);                    \
		RO12(state, b);                         \
		ADD(state, state, a, b);                \
		XOR(state, d, d, a);                    \
		RO8(state, d);                          \
		ADD(state, state, c, d);                \
		XOR(state, b, b, c);                    \
		RO7(state, b);                          \
	} while (0)

uint8_t input_state[64];
uint8_t output_state[64];
uint8_t plain_text[64];
uint8_t cipher_text[64];

uint8_t value_to_copy_out[64];

uint8_t get_key(uint8_t *k, uint8_t len)
{
	// Load key here
	memcpy(input_state + 16, k, 32);
	return 0x00;
}

uint8_t get_nonce(uint8_t *n, uint8_t len)
{
	// Load nonce here
	memcpy(input_state + 52, n, 12);
	return 0x00;
}

uint8_t get_counter(uint8_t *c, uint8_t len)
{
	// Load counter here
	memcpy(input_state + 48, c, 4);
	return 0x00;
}

uint8_t get_pt(uint8_t *pt, uint8_t len)
{
	memcpy(plain_text, pt, 64);
	return 0x00;
}

uint8_t reset(uint8_t *x, uint8_t len)
{
	// Reset key here if needed
	memset(input_state, 0, 64);
	memset(plain_text, 0, 64);
	memset(output_state, 0, 64);
	memset(cipher_text, 0, 64);

	uint32_t fixed[4] = {
		0x61707865,
		0x3320646e,
		0x79622d32,
		0x6b206574};

	memcpy(input_state, fixed, 16);

	return 0x00;
}

uint8_t read_input_state(uint8_t *x, uint8_t len)
{
	memcpy(value_to_copy_out, input_state, 64);
	simpleserial_put('r', 64, value_to_copy_out);
	return 0x00;
}

uint8_t read_plain_text(uint8_t *x, uint8_t len)
{
	memcpy(value_to_copy_out, plain_text, 64);
	simpleserial_put('r', 64, value_to_copy_out);
	return 0x00;
}

uint8_t perform_encryption(uint8_t *x, uint8_t len)
{
	/**********************************
	 * Start user-specific code here. */
	trigger_high();

	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");

	uint8_t round;
	uint8_t posn;

	/* Copy the input buffer to the output prior to the first round */
	memcpy(output_state, input_state, 64);

	/* Perform the ChaCha rounds in sets of two */
	for (round = 20; round >= 2; round -= 2)
	{
		/* Column round */
		quarterRound_indices(output_state, 0, 4, 8, 12);
		quarterRound_indices(output_state, 1, 5, 9, 13);
		quarterRound_indices(output_state, 2, 6, 10, 14);
		quarterRound_indices(output_state, 3, 7, 11, 15);

		/* Diagonal round */
		quarterRound_indices(output_state, 0, 5, 10, 15);
		quarterRound_indices(output_state, 1, 6, 11, 12);
		quarterRound_indices(output_state, 2, 7, 8, 13);
		quarterRound_indices(output_state, 3, 4, 9, 14);
	}

	// /* Add the original input to the final output */

	for (posn = 0; posn < 16; ++posn)
	{
		ADD_END(output_state, input_state, posn, posn);
	}

	/*Added by me*/
	for (posn = 0; posn < 64; ++posn)
	{
		cipher_text[posn] = output_state[posn] ^ plain_text[posn];
	}

	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	__asm("NOP\n\t");
	trigger_low();

	/* End user-specific code here. *
	********************************/
	return 0x00;
}

uint8_t read_final_state(uint8_t *x, uint8_t len)
{
	memcpy(value_to_copy_out, output_state, 64);
	simpleserial_put('r', 64, value_to_copy_out);
	return 0x00;
}

uint8_t read_ciphertext(uint8_t *x, uint8_t len)
{
	memcpy(value_to_copy_out, cipher_text, 64);
	simpleserial_put('r', 64, value_to_copy_out);
	return 0x00;
}

int main(void)
{
	platform_init();
	init_uart();
	trigger_setup();

	/* Uncomment this to get a HELLO message for debug */
	/*
	putch('h');
	putch('e');
	putch('l');
	putch('l');
	putch('o');
	putch('\n');
	*/
	reset(NULL, 0);

	simpleserial_init();
	simpleserial_addcmd('k', 32, get_key);
	simpleserial_addcmd('n', 12, get_nonce);
	simpleserial_addcmd('c', 4, get_counter);
	simpleserial_addcmd('p', 64, get_pt);
	simpleserial_addcmd('x', 0, reset);
	simpleserial_addcmd('i', 0, read_input_state);
	simpleserial_addcmd('m', 0, read_plain_text);
	simpleserial_addcmd('e', 0, perform_encryption);
	simpleserial_addcmd('f', 0, read_final_state);
	simpleserial_addcmd('o', 0, read_ciphertext);

	while (1)
		simpleserial_get();
}
