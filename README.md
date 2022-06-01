# The virtual aging brain

This repository contains the code to perform and implement the Virtual Aging Brain (VAB) as described in [1].

## Git-repo

1. Clone this repo in your laptop
   
   1. Open the terminal and run the command
   
   2. ```
      git clone https://github.com/ins-amu/virtual_aging_brain.git
      ```

2. You can always download this folder from the website

## Python Environment in Linux and MacOS

1. Open the terminal and verify to have python3 by running

```
python3
```

The version should be < 3.9 if possible

2. Run the following commands in your git folder.

```
rm -rf env
python3 -mvenv env
. env/bin/activate
pip3 install --upgrade pip
pip3 install wheel
pip3 install -r requirements.txt
pip3 install -e .
```



3. To use the environment from the Jupyter notebooks, you can either [register the kernel](https://ipython.readthedocs.io/en/stable/install/kernel_install.html), or install the `jupyter-lab` in the environment. Here we show the latter. Open a new tab on the terminal and run the following commands

```
. env/bin/activate
pip3 install jupyterlab
jupyter-lab
```

4. Copy the paste the locahost link in your favourite browser or use the deafult opening session
5. Test the tutorials (.ipynb file) by pressing Shift + Enter with your keyboard
6. Enjoy

## Python Environment in Windows

1. Create an Anaconda environment

2. Download an anaconda from the website
   
   1. https://www.anaconda.com/products/individual-d

3. Download the zip folder from the github website and unzip it. We suggest to place it in the documents folder

4. Open the Anaconda prompt

5. Move to the folder of the project by running the command (modify accordingly)

6. ```
   cd C:\Users\<Username>\Documents\eitn_school_2021_tvb_project
   ```

7. Run the following commands in the terminal

8. ```
   conda create -n env python=3.6
   conda activate env
   pip3 install wheel
   pip3 install -r requirements.txt
   pip3 install -e .
   ```
   Again, either install the kernel to global jupyter instance, or install jupyterlab in the environment.

9. Run the following command:

10. ```
    conda deactivate env
    conda activate env
    jupyter-lab
    ```

11. Copy the paste the locahost link in your favourite browser or use the deafult opening session

12. Test the tutorials (.ipynb file) by pressing Shift + Enter with your keyboard

13. Enjoy

## Python Environment in Windows Linux Subsystem

1. Windows Linux Subsystem (WSL): You can repeat the step above in the new WSL terminal if installed in windows --> See [Install WSL | Microsoft Docs](https://docs.microsoft.com/en-us/windows/wsl/install-win10#simplified-installation-for-windows-insiders)
2. Enjoy again

## License

This work is licensed under MIT license. See LICENSE.txt for the full text.

## Authors

Mario Lavanga (1), Johanna Stumme (2,3), Bahar Hazal Yalcinkaya (1), Jan Fousek (1), Christiane Jockwitz (2,3), Hiba Sheheitli (1), Nora Bittner (2,3), Meysam Hashemi (1), Spase Petkoski (1), Svenja Caspers (2,3), Viktor Jirsa (1)

(1)  Institut de Neurosciences des Systèmes (INS), Inserm, Aix-Marseille University, Marseille 13005, France 

(2) Institute of Neuroscience and Medicine (INM-1), Research Centre Jülich, Jülich, Germany

(3) Institute for Anatomy I, Medical Faculty & University Hospital Düsseldorf, Heinrich-Heine University Düsseldorf, Düsseldorf, Germany

## References

[1] M. Lavanga, J. Stumme, B. H. Yalcinkaya, J. Fousek, C. Jockwitz, H. Sheheitli, N. Bittner, M. Hashemi, S. Petkoski, S. Caspers, V. Jirsa, The virtual aging brain: a model-driven explanation for cognitive decline in older subjects (2022), bioRxiv, https://doi.org/10.1101/2022.02.17.480902 
