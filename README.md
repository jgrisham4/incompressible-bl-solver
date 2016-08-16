# incompressible-bl-solver

This code solves the incompressible boundary layer equations for constant pressure gradient flows.  It does so using one of two approaches:

 1. An explicit scheme developed by Wu (1961).
 2. Crank-Nicolson.

The options can be changed by commenting out subroutine calls
in solver.f95.  This code has only been tested on Linux machines, 
but it should work on other operating systems.  For Macintosh 
machines, make sure you have up-to-date GNU compilers installed.
For Windows users, your best bet is Cygwin (https://www.cygwin.com/).
Make sure you install GNU compilers.  The code can be compiled by
typing "make" in a terminal.  

An additional Python script is included which compares the finite difference
solution to the Blasius profile.  
