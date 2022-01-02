import configparser
from pathlib import Path

class Project:
    def __init__(self, projectDir) -> None:
        self.pdir = Path(projectDir)

    def get_project_init_file(self):
        """
        Fetch project.ini + read input
        """
        try:
            file = [f for f in self.pdir.iterdir() if f.suffix == '.ini'][0]
            config = configparser.ConfigParser()
            config.read(Path(file))
            try:
                tvFile = config['MAIN']['TreeViewerFile']
                chromFile = config['MAIN']['ChromLengths']
                gffgtfFiles = config['ADDITIONAL']['GFF_GTF'].split(" ")
            except KeyError:
                return "Could not find correct headers - project.ini malformed"
            return tvFile, chromFile, gffgtfFiles
        except IndexError:  # No project.ini file found
            msg = "No project.ini file found"
            return msg
    
    def get_gff_gtf_file(self):
        """
        Fetch gff
        """
        try:
            files = [f for f in self.pdir.iterdir() if (f.suffix == '.gff') or (f.suffix == '.gtf')]
            msg = [f"{f.name}" for f in files]
            return msg, files
        except IndexError:  # No project.ini file found
            msg = "No project.ini file found"
            return msg, False



testDir = (R"C:\Users\ajhar\OneDrive\Desktop\testProject")
testProject = Project(testDir)
testData = testProject.get_project_init_file()
testGenes = testProject.get_gff_gtf_file()

# print(testData)
# print(testGenes)

# nulltestProject = Project()
# nulltestData = nulltestProject.get_project_init_file()
# nulltestGenes = nulltestProject.get_gff_gtf_file()

# print(nulltestData)
# print(nulltestGenes)