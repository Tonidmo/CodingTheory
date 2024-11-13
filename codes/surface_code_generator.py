import numpy as np
from qecsim.models.rotatedplanar import RotatedPlanarCode

def generate_surface_code_pcm(distance : int) -> np.array:
    myCode = RotatedPlanarCode(distance, distance)
    return myCode.stabilizers

if __name__=='__main__':
    for d in [3,5,7,9,11,13,15]:
        a = generate_surface_code_pcm(d)
        # Save the matrix to a CSV file
        np.savetxt(f"pcm_matrices/distance_{d}_surface_code.csv", a, delimiter=",", fmt="%d")
        
