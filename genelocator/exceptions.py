"""Define common error types that can be checked by caller code"""


class BaseGeneLocatorException(Exception):
    DEFAULT_MESSAGE: str

    def __init__(self, message=None, *args):
        super().__init__(*args)
        self.message = message or self.DEFAULT_MESSAGE

    def __str__(self):
        return str(self.message)


class BadCoordinateException(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "The user has requested an invalid position (eg chrom for which no lookup data is available)"


class NoResultsFoundException(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "No results found"


class LookupCreateError(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "Could not generate a lookup because the raw data did not match the expected format"
