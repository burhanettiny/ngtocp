import streamlit as st

# ng cinsinden verilen miktar ile kopya sayısı hesaplama fonksiyonu
def ng_to_copy_number(ng, molar_mass):
    avogadro_number = 6.022e23  # Avogadro sayısı
    ng_in_grams = ng * 1e-9  # ng'i gram cinsine çeviriyoruz
    copy_number = (ng_in_grams * avogadro_number) / molar_mass  # Kopya sayısı hesaplama
    return copy_number

# Streamlit uygulaması
st.title("DNA/RNA Kopya Sayısı Hesaplayıcı")

# Kullanıcıdan DNA veya RNA baz dizisi girmesini isteyelim
sequence_input = st.text_area("DNA/RNA baz dizisini girin", height=200)

# Kullanıcıdan DNA/RNA türü seçmesi isteniyor
molecule_type = st.radio("DNA mı yoksa RNA mı?", ("DNA", "RNA"))

# Kullanıcıdan ng cinsinden miktar girmesi isteniyor
ng_amount = st.number_input("Miktarı girin (ng)", min_value=0.0, format="%.2f")

# Kullanıcıya özel bir molar kütle değeri girme seçeneği (1000 baz çiftli DNA için)
molar_mass = st.number_input("Molar kütleyi girin (ng/baz çifti)", min_value=0.0, value=650000.0)

# Baz dizisini kontrol etme
if sequence_input:
    # Kullanıcıdan gelen baz dizisini alıyoruz
    sequence_length = len(sequence_input)
    st.write(f"Verilen dizinin uzunluğu: {sequence_length} baz.")

    # Kopya sayısını hesaplama
    if ng_amount > 0 and molar_mass > 0:
        copy_number = ng_to_copy_number(ng_amount, molar_mass)
        st.write(f"{ng_amount} ng {molecule_type} için kopya sayısı yaklaşık {copy_number:.2e} kopyadır.")
    else:
        st.write("Lütfen geçerli bir miktar ve molar kütle girin.")
else:
    st.write("Lütfen bir DNA/RNA baz dizisi girin.")

# Hesaplama formülünü alt kısımda gösterme
st.subheader("Hesaplama Formülü:")
st.write("""
Kopya sayısı hesaplama formülü:

Kopya Sayısı = (Ng * Avogadro Sayısı) / Molar Kütle

Açıklamalar:
- Ng: Numune miktarı (ng)
- Avogadro Sayısı: 6.022 × 10²³ (Bir moldeki partikül sayısı)
- Molar Kütle: DNA/RNA'nın bir baz çifti başına kütlesi (ng/baz çifti)
""")
