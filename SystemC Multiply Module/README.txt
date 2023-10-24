
ECPS203	- 	Embedded Systems Modeling and Design
Author	-	Abyukth Kumar
Date	-	10/24/2023

Introduction to SystemC language and simulation


Brief description of the model:

We have three modules enclosed in a Top module:
	1. 	Stimulus - Decides the values to be fed to the module in test
	2. 	Module in test (Mult) - Multiplies two numbers a and b
	3. 	Monitor - Displays the multiplied output and checks if the multiplication is done correctly
The stimulus and monitor is connected to a testclk for synchronization.
All the three modules are connected by channels and send/receive data to/from channels using ports.


Issues faced while compiling the code:
	
	1.	After I integrated the code from the slides, when I compiled, the signals a and b were not recognized in top.h
	
	2.	After initializing the signals asig, bsig, fsig and clk in mon.h, when I compiled the code, I got an error saying clk.neg() - no member exists
	
	3.	The following was the output I obtained when I compiled the program after including all .cpp and .h files:
		Value of 0 * 0 	=	0     ------> Extra print 
		Value of 1 * 42	=	42
		Value of 2 * 21	=	42
		Value of 3 * 14	=	42
		Value of 6 * 7	=	42
		Value of 7 * 6	=	42
		Value of 14 * 3 =	42
		Value of 21 * 2 =	42
		Value of 42 * 1 =	42

Solutions to these issues:
	
	1.	The names were initialized in stim.h as 'A' and 'B'. Once I changed them to lowercase, it compiled successfully
	
	2.	I initialized clk as sc_in<int> instead of sc_in<bool>. Since clk has only two states, when initialized to bool, it worked
	
	3.	I had used SC_METHOD() in mon.h instead of SC_THREAD(). Once I changed that, the extra print did not occur

