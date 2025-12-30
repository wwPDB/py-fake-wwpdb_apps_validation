import argparse
import logging
#from wwpdb.apps.validation.src.pdboi.pdbdata.mmcifapiconnector import MMcifCategoryNotFoundException, MMcifAPIconnector
#from wwpdb.apps.validation.src.pdboi.pdbents.entry import Entry
import xml.etree.ElementTree as ET
from wwpdb.io.file.mmCIFUtil import mmCIFUtil


class GenerateMinimalCif:
    """
    Generate a minimal data content mmcif from an emdb_xml for the validation pipeline
    Intended to be superseded by the xml -> cif translator
    """

    def __init__(self, emdb_xml):
        """
        :param str emdb_xml: file address of the emdb_xml
        """
        self.xml_data = self.parse_xml(emdb_xml)
   
    def write_out(self, output_cif):
        """
        Writes self.xml_data into a minimal cif using self.write_minimal_cif
        :param output_cif:
        :return:

        :raises AssertionError: When xml_data not populated
        """
        assert self.xml_data
        self.write_minimal_cif(output_cif=output_cif, output_data=self.xml_data)

    def get_applied_symmetry_type(self):
        final_reconstruction = self.root.find('structure_determination_list/structure_determination/*/final_reconstruction')
        if final_reconstruction is None:
            return 'Not Provided'
        if final_reconstruction.find('applied_symmetry/point_group') is not None:
            return 'POINT'
        elif final_reconstruction.find('applied_symmetry/helical_parameters') is not None:
            return 'HELICAL'
        elif self.root.find('structure_determination_list/structure_determination/*/crystal_parameters') is not None:
            # Is 2D or 3D crystal
            if self.root.find('structure_determination_list/structure_determination/*/crystal_parameters/plane_group') is not None:
                return '2D CRYSTAL'
            elif self.root.find('structure_determination_list/structure_determination/*/crystal_parameters/space_group') is not None:
                return '3D CRYSTAL'
            else:
                logging.error("Could not find space or plane group for crystal entry")
                return 'Not Provided'
        else:
            return 'Not Provided'

    def get_authors(self):
        def switch_author_styling(author_string):
            split_author = author_string.split(' ')
            other = split_author[-1]
            last = ' '.join(split_author[:-1])
            other_dotted = '.'.join([x for x in other]) + '.'
            name = last + ', ' + other_dotted
            return name
        authors_list = self.root.find('admin/authors_list')
        authors = [switch_author_styling(author.text) for author in authors_list]
        return authors

    def find_text_or_return_empty(self, root, query):
        try:
            element = root.find(query)
            if element is not None:
                return element.text
            else:
                return ''
        except:
            return ''

    def parse_xml(self, emdb_xml):
        """
        Parses the xml file into a dictionary representation of mmcif
        :param str emdb_xml: file address of the xml
        :return: The relevant data
        :rtype dict:
        """
        logging.debug(emdb_xml)
        tree = ET.parse(emdb_xml)
        root = tree.getroot()
        self.root = root
        logging.debug(self.root)
        emdb_title = self.root.find('admin').find('title').text
        deposition_date = self.root.find('admin').find('key_dates').find('deposition').text
        
        emdb_id = self.root.attrib['emdb_id']
        
        experimental_method = 'ELECTRON MICROSCOPY'

        experiment_information = {}

        experiment_information['em_3d_reconstruction.symmetry_type'] = self.get_applied_symmetry_type()
        symmetry_type = experiment_information['em_3d_reconstruction.symmetry_type']
        if symmetry_type == 'Not Provided':
            # TODO check how this should be filled check with john about what to do here doesn't get translated
            experiment_information['pdbx_point_symmetry.schoenflies_symbol'] = 'I'
        elif symmetry_type == 'POINT':
            experiment_information['em_single_particle_entity.point_symmetry'] = self.root.findtext('structure_determination_list/structure_determination/*/final_reconstruction/applied_symmetry/point_group')
        elif symmetry_type == 'HELICAL':
            helical_parameters = self.root.find('structure_determination_list/structure_determination/*/final_reconstruction/applied_symmetry/helical_parameters')
            experiment_information['em_helical_entity.angular_rotation_per_subunit'] = self.find_text_or_return_empty(helical_parameters, 'delta_phi')
            experiment_information['em_helical_entity.axial_rise_per_subunit'] = self.find_text_or_return_empty(helical_parameters, 'delta_z')
            experiment_information['em_helical_entity.axial_symmetry'] = self.find_text_or_return_empty(helical_parameters, 'axial_symmetry')
        elif symmetry_type == '2D CRYSTAL':
            crystal_parameters = self.root.find('structure_determination_list/structure_determination/*/crystal_parameters')
            unit_cell = crystal_parameters.find('unit_cell')
            experiment_information['em_2d_crystal_entity.length_a'] = self.find_text_or_return_empty(unit_cell, 'a')
            experiment_information['em_2d_crystal_entity.length_b'] = self.find_text_or_return_empty(unit_cell, 'b')
            experiment_information['em_2d_crystal_entity.length_c'] = self.find_text_or_return_empty(unit_cell, 'c')
            experiment_information['em_2d_crystal_entity.angle_gamma'] = self.find_text_or_return_empty(unit_cell, 'gamma')
            experiment_information['em_2d_crystal_entity.space_group_name_H-M'] = self.find_text_or_return_empty(unit_cell, 'plane_group')
        elif symmetry_type == '3D CRYSTAL':
            crystal_parameters = self.root.find('structure_determination_list/structure_determination/*/crystal_parameters')
            unit_cell = crystal_parameters.find('unit_cell')
            experiment_information['em_3d_crystal_entity.length_a'] = self.find_text_or_return_empty(unit_cell, 'a')
            experiment_information['em_3d_crystal_entity.length_b'] = self.find_text_or_return_empty(unit_cell, 'b')
            experiment_information['em_3d_crystal_entity.length_c'] = self.find_text_or_return_empty(unit_cell, 'c')
            experiment_information['em_3d_crystal_entity.angle_alpha'] = self.find_text_or_return_empty(unit_cell, 'alpha')
            experiment_information['em_3d_crystal_entity.angle_beta'] = self.find_text_or_return_empty(unit_cell, 'beta')
            experiment_information['em_3d_crystal_entity.angle_gamma'] = self.find_text_or_return_empty(unit_cell, 'gamma')
            experiment_information['em_3d_crystal_entity.space_group_name'] = self.find_text_or_return_empty(unit_cell, 'space_group')

        # get images_used / number of particles
        images_used_element = self.root.find('structure_determination_list/structure_determination/*/final_reconstruction/number_images_used')
        if images_used_element is None:
            images_used_element = self.root.find('structure_determination_list/structure_determination/*/final_reconstruction/number_subtomograms_used')
        experiment_information['em_3d_reconstruction.num_particles'] = self.find_text_or_return_empty(images_used_element, '.')

        # CTF correction
        ctf_correction = self.find_text_or_return_empty(root, 'structure_determination_list/structure_determination/*/ctf_correction/details')
        experiment_information['em_3d_reconstruction.ctf_correction_method'] = ctf_correction

        # microscope information
        microscope_information = self.root.find('structure_determination_list/*/microscopy_list/*')


        experiment_information['em_imaging.microscope_model'] = self.find_text_or_return_empty(microscope_information, 'microscope')
        experiment_information['em_imaging.accelerating_voltage'] = self.find_text_or_return_empty(microscope_information, 'acceleration_voltage')
        experiment_information['em_imaging.calibrated_defocus_min'] = self.find_text_or_return_empty(microscope_information, 'calibrated_defocus_min')
        experiment_information['em_imaging.calibrated_defocus_max'] = self.find_text_or_return_empty(microscope_information, 'calibrated_defocus_max')
        experiment_information['em_imaging.nominal_defocus_min'] = self.find_text_or_return_empty(microscope_information, 'nominal_defocus_min')
        experiment_information['em_imaging.nominal_defocus_max'] = self.find_text_or_return_empty(microscope_information, 'nominal_defocus_max')
        experiment_information['em_imaging.calibrated_magnification'] = self.find_text_or_return_empty(microscope_information, 'calibrated_magnification')
        experiment_information['em_imaging.nominal_magnification'] = self.find_text_or_return_empty(microscope_information, 'nominal_magnification')


        experiment_information['em_imaging.electron_dose'] = self.find_text_or_return_empty(microscope_information, 'image_recording_list/image_recording/average_electron_dose_per_image')


        experiment_information['em_image_recording.film_or_detector_model'] = self.find_text_or_return_empty(microscope_information, 'image_recording_list/image_recording/film_or_detector_model')

        # author information
        author_list = self.get_authors()

        try:
            resolution_element = self.root.find('structure_determination_list/structure_determination/*/final_reconstruction/resolution')
            resolution = resolution_element.text
            resolution_source = resolution_element.attrib['res_type']
        except:
            resolution = ''
            resolution_source = ''
        try:
            resolution_criteria = list(self.root.iter('resolution_method'))[0].text
        except:
            resolution_criteria = ''

        experiment_information['em_3d_reconstruction.resolution'] = resolution
        experiment_information['em_3d_reconstruction.resolution_method'] = resolution_criteria

        experiment_information['em_experiment.reconstruction_method'] = self.find_text_or_return_empty(root, 'structure_determination_list/structure_determination/method')
        reconstruction_method_conversion = {
            "singleParticle": "SINGLE PARTICLE",
            "subtomogramAveraging": "SUBTOMOGRAM AVERAGING",
            "tomography": "TOMOGRAPHY",
            "electronCrystallography": "CRYSTALLOGRAPHY",
            "helical": "HELICAL"
        }
        if experiment_information['em_experiment.reconstruction_method'] and experiment_information['em_experiment.reconstruction_method'] in reconstruction_method_conversion:
            experiment_information['em_experiment.reconstruction_method'] = reconstruction_method_conversion[experiment_information['em_experiment.reconstruction_method']]

        primary_map = list(self.root.iter('map'))[0]
        contours = primary_map.iter('contour')
        primary_contours = [x for x in contours if x.attrib['primary'] == 'true']
        if len(primary_contours) > 0:
            primary_contour = primary_contours[0]
            if primary_contour.find('level') is not None:
                contour_level = primary_contour.find('level').text
            else:
                contour_level = ''
        else:
            contour_level = ''
        map_information = {
            'entry_id': emdb_id,
            'id': '1',
            'annotation_details': '',
            'format': primary_map.attrib['format'],
            'size_kb': primary_map.attrib['size_kbytes'],
            'partition': 1,
            'contour_level': contour_level,

            'axis_order_fast': primary_map.find('axis_order').find('fast').text,
            'axis_order_medium': primary_map.find('axis_order').find('medium').text,
            'axis_order_slow':  primary_map.find('axis_order').find('slow').text,

            'cell_alpha': primary_map.find('cell').find('alpha').text,
            'cell_beta': primary_map.find('cell').find('beta').text,
            'cell_gamma': primary_map.find('cell').find('gamma').text,

            'cell_a': primary_map.find('cell').find('a').text,
            'cell_b': primary_map.find('cell').find('b').text,
            'cell_c': primary_map.find('cell').find('c').text,

            'data_type': primary_map.find('data_type').text,
            'endian_type': '',

            'dimensions_col': primary_map.find('dimensions').find('col').text,
            'dimensions_row': primary_map.find('dimensions').find('row').text,
            'dimensions_sec': primary_map.find('dimensions').find('sec').text,

            'origin_col': primary_map.find('origin').find('col').text,
            'origin_row': primary_map.find('origin').find('row').text,
            'origin_sec': primary_map.find('origin').find('sec').text,

            'pixel_spacing_x': primary_map.find('pixel_spacing').find('x').text,
            'pixel_spacing_y': primary_map.find('pixel_spacing').find('y').text,
            'pixel_spacing_z': primary_map.find('pixel_spacing').find('z').text,

            'statistics_minimum': primary_map.find('statistics').find('minimum').text,
            'statistics_maximum': primary_map.find('statistics').find('maximum').text,
            'statistics_average': primary_map.find('statistics').find('average').text,
            'statistics_std': primary_map.find('statistics').find('std').text,

            'symmetry_space_group': primary_map.find('symmetry').find('space_group').text,
            'type': 'primary map',
            'file': primary_map.find('file').text
        }

        return {
            'map_information': map_information,

            'author_list': author_list,

            'emdb_title': emdb_title,
            'deposition_date': deposition_date,
            'experimental_method': experimental_method,
            'emdb_id': emdb_id,

            'experiment_information': experiment_information,
            }
        
    @staticmethod
    def write_minimal_cif(output_cif, output_data):
        """
        Writes output_data into output_cif
        :param str output_cif: File address of the desired cif
        :param dict output_data: Data to written out
        :return: None
        """
        cif = mmCIFUtil()
        cif.AddBlock(output_data['emdb_id'])
        
        logging.debug("Adding database_2")
        cif.AddCategory('database_2', ['database_id', 'database_code'])
        cif.InsertData('database_2', [['EMDB', output_data['emdb_id']]])

        logging.debug("adding em_map")
        map_items = [
            'entry_id',
            'id',
            'annotation_details',
            'contour_level',
            'format',
            'size_kb',
            'partition',
            
            'axis_order_fast',
            'axis_order_medium',
            'axis_order_slow',

            'cell_alpha',
            'cell_beta',
            'cell_gamma',

            'cell_a',
            'cell_b',
            'cell_c',

            'data_type',
            'endian_type',

            'dimensions_col',
            'dimensions_row',
            'dimensions_sec',

            'origin_col',
            'origin_row',
            'origin_sec',

            'pixel_spacing_x',
            'pixel_spacing_y',
            'pixel_spacing_z',

            'statistics_minimum',
            'statistics_maximum',
            'statistics_average',
            'statistics_std',

            'symmetry_space_group',
            'type',
            'file'

        ]
        cif.AddCategory('em_map', map_items)
        logging.debug([output_data['map_information'][x] for x in map_items])
        cif.InsertData('em_map', [[output_data['map_information'][x] for x in map_items]])

        logging.debug("Adding em_admin")
        cif.AddCategory('em_admin', ['emd_id',
                                     'current_status',
                                     'last_update',
                                     'deposition_date',
                                     'map_release_date',
                                     'details',
                                     'title'
                                     ])
        cif.InsertData('em_admin',
                       [
                           [
                               output_data['emdb_id'],
                               'REL',
                               '2000-01-01',
                               output_data['deposition_date'],
                               '2000-01-01',
                               'Nothing here',
                               output_data['emdb_title']
                               ]
                           ])

        logging.debug("Adding exptl")
        cif.AddCategory('exptl', [
            'entry_id',
            'method'
            ])
        cif.InsertData('exptl', [
            [
                'map_only',
                'ELECTRON MICROSCOPY'
                ]
            ])

        if output_data['author_list']:
            author_list = output_data['author_list']
            logging.debug("adding em_author_list")
            cif.AddCategory('em_author_list', [
                'ordinal',
                'author'
            ])
            cif.InsertData('em_author_list',
                           [[i+1, x] for i, x in enumerate(author_list)]
                           )
        else:
            logging.warning("Could not add em_author_list")

        logging.debug("Adding em_imaging")
        cif.AddCategory('em_imaging', [
            'entry_id',
            'id',
            'microscope_model',
            'accelerating_voltage',
            'calibrated_defocus_min',
            'calibrated_defocus_max',
            'nominal_defocus_min',
            'nominal_defocus_max',
            'calibrated_magnification',
            'nominal_magnification',
            'electron_dose',
        ])
        cif.InsertData('em_imaging',[
            [
                output_data['emdb_id'],
                '1',
                output_data['experiment_information']['em_imaging.microscope_model'],
                output_data['experiment_information']['em_imaging.accelerating_voltage'],
                output_data['experiment_information']['em_imaging.calibrated_defocus_min'],
                output_data['experiment_information']['em_imaging.calibrated_defocus_max'],
                output_data['experiment_information']['em_imaging.nominal_defocus_min'],
                output_data['experiment_information']['em_imaging.nominal_defocus_max'],
                output_data['experiment_information']['em_imaging.calibrated_magnification'],
                output_data['experiment_information']['em_imaging.nominal_magnification'],
                output_data['experiment_information']['em_imaging.electron_dose']
            ]
        ])
        logging.debug("Adding em_experiment")
        cif.AddCategory('em_experiment', [
            'reconstruction_method'
        ])
        cif.InsertData('em_experiment', [
            [
                output_data['experiment_information']['em_experiment.reconstruction_method']
            ]
        ])

        logging.debug("Adding em_image_recording")
        cif.AddCategory('em_image_recording', [
           'film_or_detector_model'
        ])
        cif.InsertData('em_image_recording', [
                       [
                           output_data['experiment_information']['em_image_recording.film_or_detector_model']
                       ]
        ])

        logging.debug("Adding em_3d_reconstruction")
        cif.AddCategory('em_3d_reconstruction', [
            'ctf_correction_method',
            'num_particles',
            'symmetry_type',
            'resolution',
            'resolution_method',
        ])
        cif.InsertData('em_3d_reconstruction', [
                       [
                           output_data['experiment_information']['em_3d_reconstruction.ctf_correction_method'],
                           output_data['experiment_information']['em_3d_reconstruction.num_particles'],
                           output_data['experiment_information']['em_3d_reconstruction.symmetry_type'],
                           output_data['experiment_information']['em_3d_reconstruction.resolution'],
                           output_data['experiment_information']['em_3d_reconstruction.resolution_method']
                       ]
        ])

        logging.debug("Adding symmetry data / unit cell information")
        symmetry_type = output_data['experiment_information']['em_3d_reconstruction.symmetry_type']

        if symmetry_type == "Not Provided":
            logging.debug("Adding pdbx_point_symmetry")
            cif.AddCategory('pdbx_point_symmetry',
                            [
                                'pdbx_point_symmetry'
                            ])
            cif.InsertData('pdbx_point_symmetry', [
                           [
                               output_data['experiment_information']['pdbx_point_symmetry.schoenflies_symbol']
                           ]
            ])
        elif symmetry_type == "POINT":
            logging.debug("Adding em_single_particle_entity")
            cif.AddCategory('em_single_particle_entity',
                            [
                                'point_symmetry'
                            ])
            cif.InsertData('em_single_particle_entity', [
                           [
                               output_data['experiment_information']['em_single_particle_entity.point_symmetry']
                           ]
            ])
        elif symmetry_type == 'HELICAL':
            logging.debug("Adding em_helical_entity")
            cif.AddCategory("em_helical_entity",
                            [
                                'angular_rotation_per_subunit',
                                'axial_rise_per_subunit',
                                'axial_symmetry'
                            ])
            cif.InsertData('em_helical_entity', [
                           [
                                output_data['experiment_information']['em_helical_entity.angular_rotation_per_subunit'],
                                output_data['experiment_information']['em_helical_entity.axial_rise_per_subunit'],
                                output_data['experiment_information']['em_helical_entity.axial_symmetry']
                           ]
            ])

        elif symmetry_type == '2D CRYSTAL':
            logging.debug("Adding em_2d_crystal_entity")
            cif.AddCategory("em_2d_crystal_entity",
                            [
                                'length_a',
                                'length_b',
                                'length_c',
                                'angle_gamma',
                                'space_group_name_H-M'
                            ])
            cif.InsertData("em_2d_crystal_entity", [
                           [
                               output_data['experiment_information']['em_2d_crystal_entity.length_a'],
                               output_data['experiment_information']['em_2d_crystal_entity.length_b'],
                               output_data['experiment_information']['em_2d_crystal_entity.length_c'],
                               output_data['experiment_information']['em_2d_crystal_entity.angle_gamma'],
                               output_data['experiment_information']['em_2d_crystal_entity.space_group_name_H-M']
                           ]
            ])
        elif symmetry_type == '3D CRYSTAL':
            logging.debug("Adding em_3d_crystal_entity")
            cif.AddCategory("em_3d_crystal_entity",
                            [
                                'length_a',
                                'length_b',
                                'length_c',
                                'angle_alpha',
                                'angle_beta',
                                'angle_gamma',
                                'space_group_name',
                            ])
            cif.InsertData("em_3d_crystal_entity", [
                           [
                               output_data['experiment_information']['em_3d_crystal_entity.length_a'],
                               output_data['experiment_information']['em_3d_crystal_entity.length_b'],
                               output_data['experiment_information']['em_3d_crystal_entity.length_c'],
                               output_data['experiment_information']['em_3d_crystal_entity.angle_alpha'],
                               output_data['experiment_information']['em_3d_crystal_entity.angle_beta'],
                               output_data['experiment_information']['em_3d_crystal_entity.angle_gamma'],
                               output_data['experiment_information']['em_3d_crystal_entity.space_group_name'],
                           ]
            ])

        logging.info("Writing cif to {}".format(output_cif))
        cif.WriteCif(outputFilePath=output_cif)


if __name__ == '__main__':
    logging.basicConfig(level=0)

    parser = argparse.ArgumentParser()
    
    required_args = parser.add_argument_group('required arguments (these must be specified)')
    req = required_args.add_argument
    
    req('--cif')
    req('--emdb_xml')
    
    args = parser.parse_args()
    
    GenerateMinimalCif(emdb_xml=args.emdb_xml).write_out(output_cif=args.cif)
