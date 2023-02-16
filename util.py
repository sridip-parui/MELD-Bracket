import os


class cd:
    """Context manager for changing the current working directory.
    
       Adopting from Brian M. Hunt's answer from 
       https://stackoverflow.com/questions/431684/how-do-i-change-the\
       -working-directory-in-python/24176022#24176022
    """

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
