# This is a code snapshot to the publication "Spline-Based Rotor and Stator Optimization of a Permanent Magnet Synchronous Motor"

Arxiv version: https://arxiv.org/abs/2402.19065

Printed version: will follow

# Prerequisites

In order to run the code you need
- Matlab (tested with version 2024a)
- Matlab Optimization Toolbox
- Matlab Curve Fitting Toolbox
- Matlab Global Optimization Toolbox

No other packages need to be downloaded at this point. To provide a fully working snapshot, we already include here the code snapshots for the packages

- GeoPDEs (https://rafavzqz.github.io/geopdes/)
- The NURBS toolbox (https://de.mathworks.com/matlabcentral/fileexchange/26390-nurbs-toolbox-by-d-m-spink)
- Tensor Toolbox (https://www.tensortoolbox.org/)

which we give full credit to.


# How to run the code

In the results folder, there are several scrips that can be run. These scripts reproduce the results from the publication:
- PMSM_optimization_nonlin.m - runs the optimization (takes some time)
- PMSM_eval_init.m - evaluation of the initial geometry
- PMSM_eval_opt.m - evaluation of the optimized geometry
- PMSM_std_map.m - postprocessing of the torque std for initial and optimized geometry
- PMSM_torque_profiles.m - postprocessing of the torque profiles for the initial and optimized geometry

# License

These files are free to use for research purposes. You can edit it or modify it as your needs under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

These files are distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

# Authors: 
Michael Wiesheu, Theodor Komann, Melina Merkel, Sebastian Schöps, Stefan Ulbrich, Idoia Cortes Garcia

# Support:
michael.wiesheu@tu-darmstadt.de
komann@mathematik.tu-darmstadt.de
