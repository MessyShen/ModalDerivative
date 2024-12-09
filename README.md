# LMA/MD example project

Simple implementation for LMA and MD.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

Execute
```
./build/example_bin ./input/example_bar
```
and hopefully, you will see a GUI based on igl and ImGui.

For other materials:
```
./build/example_bin ./input/example_bar --materialType NeoHookean --youngsModulus 4000 --possionRatio 0.35 --density 25.0
```

Only StVK and NeoHookean materials are implemented in this example.
By default, we use StVK with Youngs Modulus 5000, Possion Ratio 0.4 and density 20.


## Example

![Example1](./figs/Example1_LMA.png)
![Example2](./figs/Example2_MD.png)


### Notice

* If ``lma_j <= 0`` and ``lma_i >= 1``, LMA modes are visualized.
* If ``lma_j >= lma_i >= 1``, MD modes are shown.
* Since modes are not normalized, you may need to adjust the `show_scale` parameter to achieve appropriate visualization results. Always click ``Visualize Modes`` after adjustments.
* To export, click ``Export Modes`` button. ``Export Path`` needs to be a directory.