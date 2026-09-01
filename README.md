# Stack-Simulator

Everything needed for the intern project is found here.
Found in the repository are different aspects of what we need for the full project

## To run and debug:
* Use Visual Studio Code (free and available from Windows Store).
* From Visual Studio Code use File->Open Folder to open the top level directory.
* Under Project Folder open index.html.
* Change from Text Editor to Integrated Browser (upper right corner).

## To build stand alone Windows exe:
* Install Node.js from https://nodejs.org/en/download.
* From the "Project Folder" directory:
* npm init -y
* npm install --save-dev electron electron-builder
* npm start (test run)
* npm run dist (create exe)
* exe is "Project Folder/dist/stack-simulator 1.0.0.exe"
