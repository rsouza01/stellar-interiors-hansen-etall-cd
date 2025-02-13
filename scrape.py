#!/usr/bin/env python3

import requests
from bs4 import BeautifulSoup
import os
import re


start_url="http://astro.if.ufrgs.br/evol/evolve/hansen/"
base_url="http://astro.if.ufrgs.br/evol/evolve/hansen/"
FOLDER="stellar_data"

visited = set()

def download_page(url, base_url, folder):
    if url in visited:
        return
    visited.add(url)

    print(f'download_page.url: {url}')
    print(f'download_page.folder: {folder}')

    response = requests.get(url)
    if response.status_code != 200:
        return

    soup = BeautifulSoup(response.text, 'html.parser')
    page_title = soup.title.string if soup.title else 'no_title'
    page_title = page_title.replace('/', '_').replace('\\', '_')

    print(f'page_title: {page_title}')

    if not os.path.exists(folder):
        os.makedirs(folder)
    
    with open(os.path.join(folder, f"{page_title}.html"), 'w', encoding='utf-8') as file:
        file.write(response.text)

    for link in soup.find_all(re.compile("[aA]"), href=True):
        next_url = link['href']
        # if next_url in [base_url]:
        if not next_url.startswith('http'):
            next_url = base_url + next_url

        print(f'next_url: {next_url}')


        download_page(next_url, base_url, folder)

if __name__ == "__main__":
    download_page(start_url, base_url, FOLDER)