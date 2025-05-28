import logging, os

def create_logger(level=logging.INFO):
    
    """
    Create the root logger

    :return: The root logger for the program
    """

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    
    #make sure that under this name (i.e. root) only one logger is created even if function is called many times
    if not root_logger.handlers:
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setFormatter(logging.Formatter(LOG_FORMAT))
        console.setLevel(level)
        root_logger.addHandler(console)

        # Create a file handler for log messages in the output directory for the root thread
        if os.path.exists("cryptogenotyper.log"):
            os.remove("cryptogenotyper.log")
        fh = logging.FileHandler("cryptogenotyper.log", 'a', 'utf-8')
        fh.setLevel(logging.DEBUG)  # Set the file handler's level to DEBUG
        fh.setFormatter(logging.Formatter(LOG_FORMAT))
        root_logger.addHandler(fh)

    return root_logger