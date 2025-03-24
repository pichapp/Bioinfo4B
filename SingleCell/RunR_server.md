Need:

download both files to the scratch space
in the terminal, cd to the folder where the files are located
in the .sh file, change the path to the .yml file (you can right click on the .yml file to directly copy the path)
 in the .sh file, change the scrath space path in the SBATCH commands to your own, the important thing are the endings (install_seurat_%j.log/.err)(this creates an error file that tells you if something has gone wrong)
 change the email to your own.
                                            
  ![image](https://github.com/user-attachments/assets/f0821983-b1f3-41c4-b876-5a494b603fae)
                                          
                                            
After this, in the terminal write the command: sbatch pckg_dwn_seurat.sh

The download takes about 25min
you can check wether its still going with the command squeue                                                                                                        
