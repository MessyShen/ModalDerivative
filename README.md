# LMA/MD example project

Simple implementation for Linear Model Analysis (LMA) and Modal Derivative (MD).

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

This should find and build the dependencies and create a `example_bin` binary.


## Run (example)

Go back to project directory, execute
```
./build/example_bin ./input/example_bar --numModes 5
```
and hopefully, you will see a binary file exported in ``./input/example_bar/export``. ``numModes`` stands for number of LMA modes.

For other materials/number of LMA modes:
```
./build/example_bin ./input/example_bar --materialType NeoHookean --youngsModulus 4000 --possionRatio 0.35 --density 2.0 --numModes 10
```

Only StVK and NeoHookean materials are implemented in this example.
By default, we use StVK with Youngs Modulus 5000, Possion Ratio 0.4 and density 1.0.

## Customized data

* Put ``*.node`` ``*.ele`` file in one folder.
* For constained meshes, put ``*.constrainedDoFs`` in same folder.

## Notice

* We use finite difference for approximating H in Modal derivatives.
* MD modes are more suitable for constrained meshes. For unconstrained meshes, first 6 LMA modes are removed.
* Both exported LMA and MD modes are normalized around Mass-matrix by x = x / sqrt(x^T M x). You may needs an extra scaling factor to achieve appropriate poses (something like ``pose[i] = rest_shape + 0.01 * modes[i]``)

<!-- * If ``lma_j <= 0`` and ``lma_i >= 1``, LMA modes are visualized.
* If ``lma_j >= lma_i >= 1``, MD modes are shown.
* Since modes are not normalized, you may need to adjust the `show_scale` parameter to achieve appropriate visualization results. Always click ``Visualize Modes`` after adjustments.
* To export, click ``Export Modes`` button. ``Export Path`` needs to be a directory. -->