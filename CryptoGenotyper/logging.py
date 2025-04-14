import logging

def create_logger(name, level=logging.INFO):
    
    """
    Create the logger

    :return: The root logger for the program
    """

    log = logging.getLogger(name)
    LOG_FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
    
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter(LOG_FORMAT))
    log.addHandler(console)
    log.setLevel(level)

    # Create a file handler for log messages in the output directory for the root thread
    fh = logging.FileHandler("cryptogenotyper.log", 'w', 'utf-8')
    log.addHandler(fh)  

    return log