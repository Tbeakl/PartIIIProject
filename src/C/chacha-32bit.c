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

// Also replace all the memsets and memcpys with explcit sets on copys

#include "hal.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "simpleserial.h"

// ChaCha implementation taken from https://github.com/rweather/lwc-finalists/blob/5d2b22c9ff7744be429cabda0c078ea5b7b6f79e/src/individual/ASCON_masked/aead-random.c

#define leftRotate(a, bits)                           \
	(__extension__({                                  \
		uint32_t _temp = (a);                         \
		(_temp << (bits)) | (_temp >> (32 - (bits))); \
	}))

#define leftRotate7(a) (leftRotate((a), 7))
#define leftRotate8(a) (leftRotate((a), 8))
#define leftRotate12(a) (leftRotate((a), 12))
#define leftRotate16(a) (leftRotate((a), 16))

/* Perform a ChaCha quarter round operation */
#define quarterRound(a, b, c, d)              \
	do                                        \
	{                                         \
		uint32_t _b = (b);                    \
		uint32_t _a = (a) + _b;               \
		uint32_t _d = leftRotate16((d) ^ _a); \
		uint32_t _c = (c) + _d;               \
		_b = leftRotate12(_b ^ _c);           \
		_a += _b;                             \
		(d) = _d = leftRotate8(_d ^ _a);      \
		_c += _d;                             \
		(a) = _a;                             \
		(b) = leftRotate7(_b ^ _c);           \
		(c) = _c;                             \
	} while (0)

uint32_t input_state[16];
uint32_t output_state[16];
uint32_t plain_text[16];
uint32_t cipher_text[16];

uint32_t value_to_copy_out[64];

uint8_t get_key(uint8_t *k, uint8_t len)
{
	// Load key here
	memcpy(input_state + 4, k, 32);
	return 0x00;
}

uint8_t get_nonce(uint8_t *n, uint8_t len)
{
	// Load nonce here
	memcpy(input_state + 13, n, 12);
	return 0x00;
}

uint8_t get_counter(uint8_t *c, uint8_t len)
{
	// Load counter here
	memcpy(input_state + 12, c, 4);
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
	memset(input_state, 0, 16 * sizeof(uint32_t));
	memset(plain_text, 0, 16 * sizeof(uint32_t));
	memset(output_state, 0, 16 * sizeof(uint32_t));
	memset(cipher_text, 0, 16 * sizeof(uint32_t));

	input_state[0] = 0x61707865;
	input_state[1] = 0x3320646e;
	input_state[2] = 0x79622d32;
	input_state[3] = 0x6b206574;

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

	uint8_t round;
	uint8_t posn;

	/* Copy the input buffer to the output prior to the first round */
	for (posn = 0; posn < 16; ++posn)
	{
		output_state[posn] = input_state[posn];
	}

	/* Perform the ChaCha rounds in sets of two */
	for (round = 20; round >= 2; round -= 2)
	{
		/* Column round */
		quarterRound(output_state[0], output_state[4], output_state[8], output_state[12]);
		quarterRound(output_state[1], output_state[5], output_state[9], output_state[13]);
		quarterRound(output_state[2], output_state[6], output_state[10], output_state[14]);
		quarterRound(output_state[3], output_state[7], output_state[11], output_state[15]);

		/* Diagonal round */
		quarterRound(output_state[0], output_state[5], output_state[10], output_state[15]);
		quarterRound(output_state[1], output_state[6], output_state[11], output_state[12]);
		quarterRound(output_state[2], output_state[7], output_state[8], output_state[13]);
		quarterRound(output_state[3], output_state[4], output_state[9], output_state[14]);
	}

	/* Add the original input to the final output */
	for (posn = 0; posn < 16; ++posn)
	{
		output_state[posn] += input_state[posn];
	}

	/*Added by me*/
	for (posn = 0; posn < 16; ++posn)
	{
		cipher_text[posn] = output_state[posn] ^ plain_text[posn];
	}
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

