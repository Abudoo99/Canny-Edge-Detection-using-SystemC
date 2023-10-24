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

#include "mon.h"

void Mon::display()
{
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
	sc_assert(a==1);
	sc_assert(b==42);
	
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
   sc_assert(a==2);
   sc_assert(b==21);
   
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
   sc_assert(a==3);
   sc_assert(b==14);
   
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
	sc_assert(a==6);
   sc_assert(b==7);
   
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
   sc_assert(a==7);
	sc_assert(b==6);
   
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
  	sc_assert(a==14);
   sc_assert(b==3);
   
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
  	sc_assert(a==21);
   sc_assert(b==2);
		
	wait();
	cout<<"\nValue of "<<a<<" * "<<b<<"\t=\t"<<f;
   sc_assert(a==42);
   sc_assert(b==1);
}
