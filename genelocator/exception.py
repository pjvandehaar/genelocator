"""Define common error types that can be checked by caller code"""


class BaseGeneLocatorException(Exception):
    def __init__(self, message=None, *args):
        super().__init__(*args)
        self.message = message or self.DEFAULT_MESSAGE

    def __str__(self):
        return str(self.message)


# Problems creating or using gene lookups
class LookupCreateError(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "Could not generate a lookup because the raw data did not match the expected format"


class BadCoordinateException(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "The user has requested an invalid position (eg chrom for which no lookup data is available)"


class NoResultsFoundException(BaseGeneLocatorException):
    DEFAULT_MESSAGE = "No results found"


# Problems locating the cached lookup data
class UnsupportedDatasetException(BaseGeneLocatorException):
    """unable to load data (for any reason)"""
    DEFAULT_MESSAGE = "No dataset was found for your requested search"


class NoCachedDataException(UnsupportedDatasetException):
    """Unable to load data (because we didn't download a copy yet)"""
    pass


class AssetFetchError(UnsupportedDatasetException):
    DEFAULT_MESSAGE = "Could not download or generate the requested assets"
