import streamlit as st

# ng/uL cinsinden verilen miktar ile kopya sayısını hesaplama fonksiyonu
def ng_to_cp_per_ul(ng_per_ul, sequence_length, molar_mass_per_base):
    avogadro_number = 6.022e23  # Avogadro sayısı
    ng_in_grams = ng_per_ul * 1e-9  # ng'i gram cinsine çeviriyoruz
    total_molar_mass = molar_mass_per_base * sequence_length  # DNA/RNA'nın toplam molar kütlesi
    copy_number_per_ul = (ng_in_grams * avogadro_number) / total_molar_mass  # Kopya sayısı hesaplama
    return copy_number_per_ul

# Streamlit uygulaması
st.title("DNA/RNA Kopya Sayısı Hesaplayıcı")

# Kullanıcıdan DNA veya RNA türünü seçmesini iste
molecule_type = st.radio("Molekül tipi:", ("DNA", "RNA"))

# Eğer DNA seçildiyse ss veya ds seçeneğini göster
if molecule_type == "DNA":
    strand_type = st.radio("Tek zincirli (ss) mi yoksa çift zincirli (ds) mi?", ("ss", "ds"))
    molar_mass_per_base = 660 if strand_type == "ds" else 330  # dsDNA = 660 g/mol/baz, ssDNA = 330 g/mol/baz
else:
    strand_type = "ss"  # RNA genellikle tek zincirli olduğundan otomatik olarak "ss" olarak atanır.
    molar_mass_per_base = 340  # ssRNA = 340 g/mol/baz

# Kullanıcıdan baz uzunluğu girmesini iste
sequence_length = st.number_input("Baz uzunluğunu girin (baz sayısı)", min_value=1, value=1000)

# Kullanıcıdan ng/uL cinsinden miktar girmesini iste
ng_per_ul = st.number_input("Konsantrasyonu girin (ng/µL)", min_value=0.0, format="%.2f")

# Hesaplama yap
if sequence_length > 0 and ng_per_ul > 0:
    cp_per_ul = ng_to_cp_per_ul(ng_per_ul, sequence_length, molar_mass_per_base)
    st.write(f"{ng_per_ul} ng/µL {strand_type}{molecule_type} için kopya sayısı yaklaşık {cp_per_ul:.2e} kopya/µL'dir.")
else:
    st.write("Lütfen geçerli bir baz uzunluğu ve ng/µL değeri girin.")

# Hesaplama formülünü bilimsel olarak gösterme
st.subheader("Hesaplama Formülü:")
st.write("""
Kopya sayısı hesaplama formülü:

\\[
Kopya \\ Sayısı (\\text{cp/µL}) = \\frac{ (Ng/µL) \\times (6.022 × 10^{23}) }{ \\text{Molar Kütle} \\times \\text{Baz Uzunluğu} }
\\]

Açıklamalar:
- **Ng/µL**: Numune konsantrasyonu (nanogram/µL)
- **Avogadro Sayısı**: 6.022 × 10²³ (Bir moldaki molekül sayısı)
- **Molar Kütle**: DNA/RNA’nın **bir baz başına** kütlesi (g/mol)  
  - ssDNA: **330 g/mol**  
  - dsDNA: **660 g/mol**  
  - ssRNA: **340 g/mol**  
- **Baz Uzunluğu**: DNA/RNA dizisinin uzunluğu (baz sayısı)
""")
