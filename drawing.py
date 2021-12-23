import numpy as np
import matplotlib.pyplot as plt

from utils import create_spheres, draw_sphere, plot
from classes import Sphere

params = {
    'n': 5,
    'g': 9.8,
    'radius_min': 1,
    'radius_max': 5,
    'x_min': 0,
    'x_max': 20,
    'y_min': 0,
    'y_max': 20,
    'z_min': 0,
    'z_max': 20
}

arr = [Sphere(center_coords=[10, 10, 10.1], radius=5, mass=10, g=9), Sphere(center_coords=[10, 10, 0], radius=5, mass=1, g=9)]

# plot(arr)
# plt.show()

print(arr[0].collides_with(arr[1]))
print(arr[0].distance(arr[1]))

