
import  setuptools

def readme():
    with open('README.md') as f:
        return f.read()

setuptools.setup(name='SWIT',
      version='1.0',
      description='Seismic Waveform Inversion Toolbox-1.0',
      author='Haipeng Li',
      author_email='haipengl@mail.ustc.edu.cn',
      license='GNU General Public License v3.0',
      packages=setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",],
      install_requires=['obspy',
                        'numpy',
                        'scipy',
                        'multiprocess',
                        'matplotlib', 
                        'PySimpleGUI',
                        'psutil',
                        'Pillow'],
      python_requires='>=3.7',
      zip_safe=False,
      )