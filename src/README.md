##Instructions

Cleanup:
```bash
make clean
```

Build:
```bash
make
```

Create directories for output:
```bash
make prepare
```
Run simulation:
```bash
./run
```

Create an image output of a timestep:
```bash
python image_create.py output/100.dat 5.0
```
(creates an image of the 100th timestep with 5x image brightness adjustment)




