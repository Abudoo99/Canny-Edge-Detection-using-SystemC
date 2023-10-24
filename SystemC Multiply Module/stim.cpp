/*
 * ECPS203 - Assignment 3
 *
 * Introduction to SystemC language and simulation
 *    
 * Source code adapted from "SystemC Training: The Definitive Guide to SystemC"
 * by Doulos Ltd., 2015
 *
 * Author: Abyukth Kumar
 *
*/

#include "stim.h"

void Stim::stimulus()
{
	wait();
	a = 1;
	b = 42;
	wait();
	a = 2;
	b = 21;
	wait();
	a = 3;
	b = 14;
	wait();
	a = 6;
	b = 7;
	wait();
	a = 7;
	b = 6;
	wait();
	a = 14;
	b = 3;
	wait();
	a = 21;
	b = 2;
	wait();
	a = 42;
	b = 1;
	wait();
	sc_stop();
}
