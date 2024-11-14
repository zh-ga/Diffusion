from setuptools import setup
 
setup(
    name='diffusion_sicm',
    version='0.2',
    description='Mulitlayer Chemical Vapor Deposition Diffusion Analytical solution',
    packages=['diffusion'],
    install_requires=[
        "numpy",
        "scipy",
        "pyyaml",
    ],
    python_requires='>=3.10',
)