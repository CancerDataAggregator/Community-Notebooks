import sevenbridges as sbg
import pandas as pd
from multiprocessing import Pool
import json
import time
import os
import math
from IPython.display import clear_output

pd.set_option('max_columns', None)

CGC_ENDPOINT = 'https://cgc-api.sbgenomics.com/v2'

CHUNK_LIMIT = 100

CPU = os.cpu_count()

INTEGER_FIELDS = [
    'Patient_days_to_birth',
    'Diagnosis_age_at_diagnosis',
    'Specimen_age_at_collection',
    'File_byte_size',
    'age_at_diagnosis'
]

CGC_MAPPING = {
        'sample_id': 'Specimen_id',
        'investigation': 'ResearchSubject_associated_project',
        'case_id': 'ResearchSubject_id',
        'primary_site': 'ResearchSubject_primary_disease_site',
        'disease_type': 'ResearchSubject_primary_disease_type',
        'gender': 'Patient_sex',
        'age_at_diagnosis': 'Diagnosis_age_at_diagnosis',
        'race': 'Patient_race',
        'ethnicity': 'Patient_ethnicity',
        'sample_type': 'Specimen_source_material_type'
    }


MAPPING_DICT = {
    'id': 'File_id',
    'identifier': 'File_identifier',
    'label': 'File_label',
    'data_category': 'File_data_category',
    'data_type': 'File_data_type',
    'file_format': 'File_format',
    'associated_project': 'File_associated_project',
    'drs_uri': 'File_drs_uri',
    'byte_size': 'File_byte_size',
    'checksum': 'File_checksum',
    'data_modality': 'File_data_modality',
    'imaging_modality': 'File_imaging_modality',
    'dbgap_accession_number': 'File_dbgap_accession_number',
    'subject_id': 'Patient_id',
    'sex': 'Patient_sex',
    'race': 'Patient_race',
    'ethnicity': 'Patient_ethnicity',
    'days_to_birth': 'Patient_days_to_birth',
    'researchsubject_id': 'ResearchSubject_id',
    'member_of_research_project': 'ResearchSubject_associated_project',
    'primary_diagnosis_condition': 'ResearchSubject_primary_disease_type',
    'primary_diagnosis_site': 'ResearchSubject_primary_disease_site',
    'researchsubject_specimen_id': 'Specimen_id',
    'associated_project_specimen': 'Specimen_associated_project',
    'days_to_collection': 'Specimen_days_to_collection',
    'primary_disease_type': 'Specimen_primary_disease_type',
    'anatomical_site': 'Specimen_anatomical_site',
    'source_material_type': 'Specimen_source_material_type',
    'specimen_type': 'Specimen_type',
    'derived_from_specimen': 'Specimen_derived_from_specimen',
    'primary_diagnosis': 'Disease_primary_diagnosis',
    'age_at_diagnosis': 'Disease_age_at_diagnosis',
    'morphology': 'Disease_morphology',
    'stage': 'Disease_stage',
    'grade': 'Disease_grade',
    'method_of_diagnosis': 'Disease_method_of_diagnosis',
    'treatment_type': 'Treatment_treatment_type',
    'treatment_outcome': 'Treatment_treatment_outcome',
    'days_to_treatment_start': 'Treatment_days_to_treatment_start',
    'days_to_treatment_end': 'Treatment_days_to_treatment_end',
    'therapeutic_agent': 'Treatment_therapeutic_agent',
    'treatment_anatomic_site': 'Treatment_treatment_anatomic_site',
    'treatment_effect': 'Treatment_treatment_effect',
    'treatment_end_reason': 'Treatment_treatment_end_reason',
    'number_of_cycles': 'Treatment_number_of_cycles',
}


def iter_pages(result):
    """Iterates over pages of a query result to return results as a single Dataframe.
    
    Inputs:
        results: object of class Result. Obtained by running a query.
    Returns:
        Pandas Dataframe containing query results.
    """
    df = result.to_dataframe()
    while result.has_next_page:
        result = result.next_page()
        df = pd.concat([df, result.to_dataframe()])
    return df

def log_state(import_jobs):
    """Logs state of the bulk import process.
    """
    state_list = [job.state for job in import_jobs]
    s = state_list.count('SUBMITTED')
    r = state_list.count('RUNNING')
    f = state_list.count('FINISHED')
    t = len(import_jobs)
    print(f"\nSUBMITTED: {s}")
    print(f"RUNNING: {r}")
    print(f"FINISHED: {f}/{t}")


def process_and_upload(df, token, endpoint=CGC_ENDPOINT, import_project=None, tags=None, conflict_resolution='SKIP'):
    """Imports CDA files to a CGC project from a given pandas Dataframe.
    """
    
    # Prepare dataframe for processing
    metadata_fields = list(MAPPING_DICT.values())
    metadata_fields.remove('File_identifier')
    
    df.rename(columns=MAPPING_DICT, inplace=True)

    # Preparing upload chunks:
    files_to_upload = {}
    for i, row in df.iterrows():
        if type(row.File_id) != str:
            continue
        if row.File_id in files_to_upload:
            for field in df.columns:
                try:
                    if type(row[field]) == float:
                        if math.isnan(row[field]):
                            if type(files_to_upload[row.File_id]['metadata'][field]) == float:
                                if not math.isnan(files_to_upload[row.File_id]['metadata'][field]):
                                    files_to_upload[row.File_id]['metadata'][field] = 'MULTIPLE'
                            else:
                                files_to_upload[row.File_id]['metadata'][field] = 'MULTIPLE'
                        else:
                            if files_to_upload[row.File_id]['metadata'][field] != row[field]:
                                files_to_upload[row.File_id]['metadata'][field] = 'MULTIPLE'
                    else:
                        if files_to_upload[row.File_id]['metadata'][field] != row[field]:
                            files_to_upload[row.File_id]['metadata'][field] = 'MULTIPLE'
                except:
                    print(field)
                    print('-----')
                    print(row)
                    print('-----')
                    print(files_to_upload[row.File_id])
                    raise(KeyError)
        else:
            files_to_upload[row.File_id] = {
                'name': row.File_label,
                'drs_uri': row.File_drs_uri,
                'project': import_project,
                'metadata': {field: row[field] for field in df.columns}
            }
            
    # Fix nan values
    for file in files_to_upload:
        for field in df.columns:
            if type(files_to_upload[file]['metadata'][field]) == float:
                if math.isnan(files_to_upload[file]['metadata'][field]):
                    files_to_upload[file]['metadata'][field] = None
        
            
    # Create chunks:
    cda_chunk = []
    cda_import = []
    for file_id in files_to_upload:
        cda_chunk.append(files_to_upload[file_id])
        if len(cda_chunk) >= CHUNK_LIMIT:
            cda_import.append(cda_chunk)
            cda_chunk = []
    
    if cda_chunk:
        cda_import.append(cda_chunk)
    
    # Connect to SBG API
    api = sbg.Api(url=endpoint, token=token)
    import_jobs = []

    expected_chunk_number = len(cda_import)

    # Initiating import jobs
    print(f"Submitting files in chunks of {CHUNK_LIMIT}.\n")
    for n, chunk in enumerate(cda_import):
        print(f"Importing chunk {n + 1}/{expected_chunk_number}")
        try:
            import_jobs.append(
                api.drs_imports.bulk_submit(imports=chunk, tags=tags, conflict_resolution=conflict_resolution))
        except sbg.errors.SbgError as e:
            time.sleep(10)
            for job in import_jobs:
                job.reload()
            clear_output(wait=True)
            print(f"Error occurred when submitting chunk {n}.")
            print(f"Reason:\n{e.message}\n{e.more_info}")
            print(f"Current status:")
            log_state(import_jobs)
            # print(traceback.format_exc())
            return chunk

    # Logging import progress
    for job in import_jobs:
        job.reload()
    state_list = [job.state for job in import_jobs]
    while state_list.count('FINISHED') != len(import_jobs):
        log_state(import_jobs)
        time.sleep(10)
        for job in import_jobs:
            job.reload()
        state_list = [job.state for job in import_jobs]
    log_state(import_jobs)
    print("\nImport completed!")
    return import_jobs