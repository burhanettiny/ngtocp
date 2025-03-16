import streamlit as st

# ng cinsinden verilen miktar ile kopya sayısı hesaplama fonksiyonu
def ng_to_copy_number(ng_per_ul, molar_mass, sequence_length):
    avogadro_number = 6.022e23  # Avogadro sayısı
    ng_in_grams = ng_per_ul * 1e-9  # ng'i gram cinsine çeviriyoruz
    total_mass = molar_mass * sequence_length  # Toplam molar kütleyi hesaplıyoruz
    copy_number_per_ul = (ng_in_grams * avogadro_number) / total_mass  # Kopya sayısı hesaplama
    return copy_number_per_ul

# Streamlit uygulaması
st.title("DNA/RNA Kopya Sayısı Hesaplayıcı")

# Kullanıcıdan DNA veya RNA baz dizisi girmesini isteyelim
sequence_input = st.text_area("DNA/RNA baz dizisini girin", height=200)

# Kullanıcıdan DNA/RNA türü seçmesi isteniyor
molecule_type = st.radio("DNA mı yoksa RNA mı?", ("DNA", "RNA"))

# Kullanıcıdan baz uzunluğu girmesi isteniyor
sequence_length = st.number_input("Baz dizisinin uzunluğunu girin (baz sayısı)", min_value=1, value=1000)

# Kullanıcıdan ng/uL cinsinden miktar girmesi isteniyor
ng_per_ul = st.number_input("Miktarı girin (ng/µL)", min_value=0.0, format="%.2f")

# Kullanıcıya özel bir molar kütle değeri girme seçeneği (1000 baz çiftli DNA için)
molar_mass = st.number_input("Molar kütleyi girin (ng/baz çifti)", min_value=0.0, value=650000.0)

# Baz dizisini kontrol etme
if sequence_input:
    st.write(f"Verilen dizinin uzunluğu: {sequence_length} baz.")

    # Kopya sayısını hesaplama
    if ng_per_ul > 0 and molar_mass > 0 and sequence_length > 0:
        copy_number_per_ul = ng_to_copy_number(ng_per_ul, molar_mass, sequence_length)
        st.write(f"{ng_per_ul} ng/µL {molecule_type} için kopya sayısı yaklaşık {copy_number_per_ul:.2e} kopyadır.")
    else:
        st.write("Lütfen geçerli bir miktar, molar kütle ve baz uzunluğu girin.")
else:
    st.write("Lütfen bir DNA/RNA baz dizisi girin.")

# Hesaplama formülünü alt kısımda bilimsel olarak gösterme
st.subheader("Hesaplama Formülü:")
st.write("""
Kopya sayısı hesaplama formülü:

Kopya Sayısı per µL = (Ng/µL * Avogadro Sayısı) / (Molar Kütle * Baz Uzunluğu)

Açıklamalar:
- Ng/µL: Numune miktarı (ng/µL)
- Avogadro Sayısı: 6.022 × 10²³ (Bir moldeki partikül sayısı)
- Molar Kütle: DNA/RNA'nın bir baz çifti başına kütlesi (ng/baz çifti)
- Baz Uzunluğu: DNA/RNA dizisinin uzunluğu (baz sayısı)
""")
