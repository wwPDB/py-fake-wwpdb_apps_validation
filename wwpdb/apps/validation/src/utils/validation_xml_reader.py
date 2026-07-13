class ValidationXMLReader:
    def __init__(self, xml_file=None, xml_string=None):
        self.failed_programs = []

    def get_failed_programs(self):
        return self.failed_programs

