#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install Faker')

from faker import Faker

fake = Faker()

def generate_synthetic_sequence(length):
    """
    Generate a synthetic protein sequence of a specified length using the Faker library.

    Parameters:
    - length: The length of the synthetic sequence (default is 50).

    Returns:
    A synthetic protein sequence as a string.
    """
    synthetic_sequence = ''.join([fake.random_element("ACDEFGHIKLMNPQRSTVWY") for _ in range(length)])
    return synthetic_sequence

def generate_synthetic_dataset(num_sequences, sequence_length):
    """
    Generate a synthetic dataset of protein sequences.

    Parameters:
    - num_sequences: The number of synthetic sequences in the dataset (default is 100).
    - sequence_length: The length of each synthetic sequence (default is 50).

    Returns:
    A list containing synthetic protein sequences.
    """
    synthetic_dataset = [generate_synthetic_sequence(sequence_length) for _ in range(num_sequences)]
    return synthetic_dataset

# Example usage:
num_sequences = 10
sequence_length = 10

synthetic_dataset = generate_synthetic_dataset(num_sequences, sequence_length)

# Print the synthetic dataset
for i, sequence in enumerate(synthetic_dataset):
    print(f'Sequence {i+1}: {sequence}')


# In[ ]:




