var SETTINGS = {
    TEST_MODE: false,
    DEBUG_MODE: false,
    USE_ECHO_SERVICE: false,
    USE_FLICKR: false,
    MAX_HTML_SIZE: 1000, //Default size = 1 GB / 1000 MB
    USE_ACROBAT: false,
    ANALYTICS: false,
    SUPPORTED_VERSION: 41,
	READER_VER: 13,		//Right now we are chosing Acrobat Reader Shim's Version as 13, if we chnage this constant, please ensure we also update the code in Native Host
	IS_READER: false,
	ERP_READER_VER: 12,
	IS_ERP_READER: false,
	IS_ACROBAT: true,
    ANALYTICS_OPT_IN: true
};
SETTINGS.USE_ACROBAT = true; 
SETTINGS.ANALYTICS = true; 
