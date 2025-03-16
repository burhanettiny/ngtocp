import streamlit as st
from Bio import Entrez
from Bio import SeqIO

# Genom verisini çekme fonksiyonu
def get_genome_length(organism_name):
    Entrez.email = "your_email@example.com"  # NCBI için geçerli bir e-posta adresi girin
    
    # Organizmaların arasından uygun olanı bulmak için sorgu yapıyoruz
    search_handle = Entrez.esearch(db="nucleotide", term=organism_name, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    # Eğer sonuç bulunursa, genoma ait veriyi çekelim
    if search_results["Count"] != "0":
        # Genom verisini çekmek için, arama sonucundaki ilk kayıtla ilgili detayları alıyoruz
        seq_id = search_results["IdList"][0]
        fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        genome_record = SeqIO.read(fetch_handle, "genbank")  # GenBank formatında okuma
        fetch_handle.close()
        
        # Genom uzunluğunu almak için seq.length özelliğini kullanıyoruz
        genome_length = genome_record.seq.length  # Genom uzunluğunu doğrudan alıyoruz
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
st.title("DNA/RNA Kopya Sayısı Hesaplayıcı")

# Kullanıcıdan organizma seçmesi isteniyor
organisms = ["Homo sapiens", "Mus musculus", "Escherichia coli", "Arabidopsis thaliana"]
selected_organism = st.selectbox("Organizma Seçin", organisms)

# Kullanıcıdan DNA/RNA türü seçmesi isteniyor
molecule_type = st.radio("DNA mı yoksa RNA mı?", ("DNA", "RNA"))

# Kullanıcıdan ng cinsinden miktar girmesi isteniyor
ng_amount = st.number_input("Miktarı girin (ng)", min_value=0.0, format="%.2f")

# Kullanıcıya genom uzunluğuna ilişkin bir sayı girme seçeneği ekleniyor
user_genome_length = st.number_input("Eğer özel bir genom uzunluğu girecekseniz, buraya yazın (Baz Sayısı)", min_value=0)

# Genom uzunluğuna bağlı olarak molar kütle hesaplama
if selected_organism == "Homo sapiens":
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için
elif selected_organism == "Mus musculus":
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için
elif selected_organism == "Escherichia coli":
    molar_mass = 660000  # Ortalama 1000 baz çiftli DNA için
else:
    molar_mass = 650000  # Ortalama 1000 baz çiftli DNA için

# Eğer kullanıcı genom uzunluğunu girmemişse, o zaman NCBI'den çekilecek
if user_genome_length > 0:
    genome_length = user_genome_length
    st.write(f"Kullanıcı tarafından girilen genom uzunluğu: {genome_length} baz.")
else:
    # Seçilen organizmanın genom bilgisi çekilecek
    if selected_organism:
        try:
            genome_length = get_genome_length(selected_organism)
            if genome_length:
                st.write(f"Seçilen organizmanın genom uzunluğu: {genome_length} baz.")
            else:
                st.write(f"{selected_organism} için genom bilgisi çekilemedi.")
        except:
            st.write("Genom bilgisi çekilemedi, lütfen tekrar deneyin.")

# Kopya sayısını hesaplama
if ng_amount > 0 and genome_length > 0:
    copy_number = ng_to_copy_number(ng_amount, molar_mass)
    st.write(f"{ng_amount} ng {selected_organism} {molecule_type} için kopya sayısı yaklaşık {copy_number:.2e} kopyadır.")
else:
    st.write("Lütfen geçerli bir miktar girin.")

# Seçilen organizmanın genom büyüklüğünü alt kısımda gösterme
if selected_organism:
    st.subheader(f"Seçilen {selected_organism} organizmasının genom büyüklüğü:")
    if user_genome_length > 0:
        st.write(f"**Girilen genom uzunluğu:** {user_genome_length} baz.")
    elif genome_length:
        st.write(f"**Genom uzunluğu NCBI'den alındı:** {genome_length} baz.")
    else:
        st.write("Genom bilgisi bulunamadı.")
