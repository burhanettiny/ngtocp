import streamlit as st
from Bio import Entrez

# Genom verisini çekme fonksiyonu
def get_genome_length(organism_name):
    Entrez.email = "mailtoburhanettin@gmail.com"  # NCBI için geçerli bir e-posta adresi girin
    
    # Organizmaların arasından uygun olanı bulmak için sorgu yapıyoruz
    search_handle = Entrez.esearch(db="nucleotide", term=organism_name, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    # Eğer sonuç bulunursa, genoma ait veriyi çekelim
    if search_results["Count"] != "0":
        # Genom verisini çekmek için, arama sonucundaki ilk kayıtla ilgili detayları alıyoruz
        seq_id = search_results["IdList"][0]
        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        genome_record = fetch_handle.read()
        fetch_handle.close()
        
        # Genomun uzunluğunu bulmak için (GB formatında) dosya okuma işlemi yapıyoruz
        genome_length = genome_record.count('N')  # 'N' olanları sayabiliriz (bazen 'N' bulunabilir)
        
        return genome_length
    else:
        print(f"{organism_name} için genom verisi bulunamadı.")
        return None

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
        genome_length = get_genome_length(selected_organism)
        if genome_length:
            st.write(f"Seçilen organizmanın genom uzunluğu: {genome_length} baz çifti.")
        else:
            st.write(f"{selected_organism} için genom bilgisi çekilemedi.")
    except:
        st.write("Genom bilgisi çekilemedi, lütfen tekrar deneyin.")

# Kullanıcıdan DNA/RNA türü seçmesi isteniyor
molecule_type = st.radio("DNA mı yoksa RNA mı?", ("DNA", "RNA"))

# Kullanıcıdan ng cinsinden miktar girmesi isteniyor
ng_amount = st.number_input("Miktarı girin (ng)", min_value=0.0, format="%.2f")

# Kullanıcıya istediği uzunlukta baz dizisi girme seçeneği ekleniyor
custom_sequence = st.text_area("Eğer özel bir baz dizisi kullanmak isterseniz, buraya yazın:")

# Eğer kullanıcı özel bir baz dizisi girdiyse, uzunluğunu hesaplayalım
if custom_sequence:
    custom_sequence_length = len(custom_sequence)
    st.write(f"Girdiğiniz özel baz dizisinin uzunluğu: {custom_sequence_length} baz.")
    genome_length = custom_sequence_length  # Kullanıcının verdiği uzunluğu kullan

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
