import shutil
import os

source = "/home/stardust/Documents/Fusion_TWAS/WEIGHTSold"
destination = "/home/stardust/Documents/Fusion_TWAS/WEIGHTS"

dirs = os.listdir(source)
# open each dir in dirs, move all files in each dir to destination
for dir in dirs:
    dir_path = os.path.join(source, dir)
    if os.path.isdir(dir_path):
        files = os.listdir(dir_path)
        for file in files:
            file_path = os.path.join(dir_path, file)
            if os.path.isfile(file_path):
                shutil.move(file_path, destination)