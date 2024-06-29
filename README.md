# bifurcation-analyser
Uses CAPD (Computer Assisted Proofs in Dynamics) to validate the number of zeroes of a map over a compact region in parameter space.

## Setup:
Create a `{path}/components` directory, where path is the directory that was cloned to, and compile using `make`.

## Usage:
Remove any pre-existing files in `{path}/components`:
```bash
rm {path}/components/*
```
After building with CMake, run
```bash
./bifurcation_order <param1_lb param1_ub> [param2_lb param2_ub...]
```
The results are sent to standard output.

If there are exactly two parameters, the components directory will be populated with files containing regular, special, and extra special boxes. Use
```bash
python3 rectangles.py
```
to display these boxes.
- Regular boxes (displayed in shades of blue per connected component) are boxes that were verified by `bisection` to contain no bifurcations. See the output of `newton` for the number of zeroes in each connected component.
- Verified boxes (displayed in red) were not able to be eliminated by `bisection`, but were later verified by `bifurcation_order` to have at most `max_derivative` zeroes.
- Special boxes (displayed in yellow) were not able to be verified to have less than `max_derivative` zeroes.
