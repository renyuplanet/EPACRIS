# EPACRIS
Version 1.0

EPACRIS (ExoPlanet Atmospheric Chemistry & Radiative Interaction Simulator) is a tool to compute the pressure-temperature profile and chemical composition of exoplanet atmospheres, based on fundamental principles and user-specified initial and boundary conditions.

EPACRIS has the capability of simulating massive atmospheres of gaseous planets, thin atmospheres of terrestrial planets and exoplanets (e.g., Earth, Mars, Venus, and rocku exoplanets), and intermediate atmospheres of sub-Neptunes, the most populous type of planets in the Galaxy.

EPACRIS is coded in C.

The current version includes the climate module of EPACRIS, i.e., the computation of pressure-temperature profiles in radiative-convective equilibrium.

Future release will include the chemistry module.

## Authors
* [Renyu Hu](https://renyuplanet.github.io/) (Jet Propulsion Laboratory, California Institute of Technology)
* Markus Scheucher (Jet Propulsion Laboratory, California Institute of Technology)

## Usage
* To compile: gcc epacris.c and to run ./a.out
* Parameter files are typically planet_.h. Each main_.c will include one planet_.h file.

## Acknowledgement
The research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004).

## License
Copyright © 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

Licensed under the Apache License, Version 2.0 (the "Licence");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## BUGS!!!
For any issues and bugs please send an e-mail at [renyu.hu@jpl.nasa.gov](renyu.hu@jpl.nasa.gov), or submit an issue through the Github system.
