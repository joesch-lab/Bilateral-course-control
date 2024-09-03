import tkinter as tk
from tkinter import filedialog
import os

def file_chooser():
    application_window = tk.Tk()

    # Build a list of tuples for each file type the file dialog should display
    my_filetypes = [('all files', '.*'), ('text files', '.txt')]

    # Ask the user to select a folder.
    answer = filedialog.askdirectory(parent=application_window,
                                     initialdir=os.getcwd(),
                                     title="Please select a folder:")
    print(answer)
    return answer

if __name__ == "__main__":
    print(file_chooser())