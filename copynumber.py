import streamlit as st
from Bio import Entrez

# NCBI'den genom bilgisi çekme fonksiyonu
def fetch_genome_info(organism):
    Entrez.email = "mailtoburhanettin@gmail.com"  # Biopython için email adresinizi girin
    search_handle = Entrez.esearch(db="nucleotide", term=organism + "[Orgn]", retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    # İlk sonucun accession numarasını al
    accession_number = search_results["IdList"][0]
    
    # Accession numarasına ait detayları çek
    fetch_handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = Entrez.read(fetch_handle)
    fetch_handle.close()
    
    # Genom uzunluğunu (baz çiftleri olarak) çek
    genome_length = int(record[0]['GBSeq_length'])
    
    return genome_length

# ng cinsinden verilen miktar ile kopya sayısı hesaplama fonksiyonu
def ng_to_copy_number(ng, molar_mass):
    avogadro_number = 6.022e23
    ng_in_grams = ng * 1e-9
    copy_number = (ng_in_grams * avogadro_number) / molar_mass
    return copy_number

# Streamlit uygulaması
st.title("DNA/RNA Kopya Sayısı Hesaplayıcı ve Genom Bilgisi")

# Kullanıcıdan organizma seçmesi isteniyor
organisms = ["Homo sapiens", "Mus musculus", "Escherichia coli", "Arabidopsis thaliana"]
selected_organism = st.selectbox("Organizma Seçin", organisms)

# Seçilen organizmanın genom bilgisini çekme
if selected_organism:
    try:
        genome_length = fetch_genome_info(selected_organism)
        st.write(f"Seçilen organizmanın genom uzunluğu: {genome_length} baz çifti.")
    except:
        st.write("Genom bilgisi çekilemedi, lütfen tekrar deneyin.")

# Kullanıcıdan DNA/RNA türü seçmesi isteniyor
molecule_type = st.radio("DNA mı yoksa RNA mı?", ("DNA", "RNA"))

# Kullanıcıdan ng cinsinden miktar girmesi isteniyor
ng_amount = st.number_input("Miktarı girin (ng)", min_value=0.0, format="%.2f")

# Genom uzunluğuna bağlı olarak molar kütle hesaplama
if selected_organism == "Homo sapiens":
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için
elif selected_organism == "Mus musculus":
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için
elif selected_organism == "Escherichia coli":
    molar_mass = 660000  # Ortalama 1000 baz çiftli DNA için
else:
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için

# Kopya sayısını hesaplama
if ng_amount > 0 and genome_length > 0:
    copy_number = ng_to_copy_number(ng_amount, molar_mass)
    st.write(f"{ng_amount} ng {selected_organism} {molecule_type} için kopya sayısı yaklaşık {copy_number:.2e} kopyadır.")
else:
    st.write("Lütfen geçerli bir miktar girin.")
