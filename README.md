# Academic Scarping

The purpose of this project is to enable easy scraping of information from academic journals.

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)


## Table of Contents

- [About](#about)
- [Setup](#setup)
- [Usage](#usage)


## About

### (1) Getting Author Emails

The first script [get_author_emails.py](./get_author_emails.py) enables you to find the emails of authors who have published in an area of your interest.

You can run this scripe as specified [below](#get_author_emailspy).

### (2) Planned future work
The next planned script is to enable the download of academic paper content and abstract. The main anticipated utility of this is for NLP, including the use of large language models.

All other script suggestions are warmly invited. (Please submit an Issue or Pull Request.)


## Setup
### Requirements
Requirements are specified in [requirements.txt](./requirements.txt) and can be installed with the following command:

```
pip install -r requirements.txt
```

### Email
To prevent the PubMed API from restricting your number of requests, you can add an email to the [config.py](./config.py) file. This will be used by the BioPython API when sending your requests.



## Usage

### [get_author_emails.py](./get_author_emails.py)
To generate a .csv file of author emails in the [output folder](./output/), run the script as follows:

```
python3 get_author_emails.py "search term" --affiliation "optional affiliation term"
```

For example:
```
python3 get_author_emails.py "proteomics nutrition" --affiliation "GSK"
```

will return all papers related to both the keywords "proteomics" and "nutrition" and where atleast one author on the paper has an affiliation with GSK.
