class Error(Exception):
    """Base class for other exceptions"""
    pass
class ConfigNotFound(Error):
    """Raise when config file is required but not present"""
    print(f"ERROR:FileNotFound: Please provide configuration file using '-c' argument or check file path")
    exit(1)