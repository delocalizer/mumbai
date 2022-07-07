from setuptools import setup

setup (
    entry_points={
        'console_scripts': [
            'mumbai = mumbai.mumbai:main',
            'mumbai_db = mumbai.mumbai_db:main',
         ]
    }
)
