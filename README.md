# psto

This project wraps a 2D Finite Element code within a heuristic optimization algorithm to solve the volume-constrained material placement problem. Written in support of my master's thesis.

## Getting Started

The `psto.m` script is the highest-level code and takes user inputs. The `particle.m` and `population.m` files define classes that call functions held in the *fe_functions*, *opt_functions*, and *view_functions* subdirectories.

## License

Copyright (C) 2020  Matt Ireland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/).