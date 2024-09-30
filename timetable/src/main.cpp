/******************************************************************************
* Copyright (C) Siarhei Uzunbajakau, 2023.
*
* This program is free software. You can use, modify, and redistribute it under
* the terms of the GNU Lesser General Public License as published by the Free
* Software Foundation, either version 3 or (at your option) any later version.
* This program is distributed without any warranty.
*
* Refer to COPYING.LESSER for more details.
******************************************************************************/

// Tests the time tables.

#define BOOST_ALLOW_DEPRECATED_HEADERS

#include <deal.II/base/timer.h>

#define TMR(__name) \
	TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

void test()
{
	double d = 0.0;
	for (unsigned int i = 0; i < 200000; i++)
		d = i*i + d;
}

int main()
{

	TimerOutput::OutputFrequency tf = TimerOutput::summary;

	TimerOutput timer(
		std::cout,
		tf,
		TimerOutput::cpu_and_wall_times_grouped);

	{TMR("Make mesh"); test();}
	{TMR("Fill Dirichlet stack"); test();}
	{TMR("Setup"); test();}
	{TMR("Assemble"); test();}
	{TMR("Solve"); test();}

return 0;
}

