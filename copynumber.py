import streamlit as st
import numpy as np

# ng/uL cinsinden verilen miktar ile kopya sayısını hesaplama fonksiyonu
def ng_to_cp_per_ul(ng_per_ul, sequence_length, molar_mass_per_base):
    avogadro_number = 6.022e23  # Avogadro sayısı
    ng_in_grams = ng_per_ul * 1e-9  # ng'yi gram cinsine çeviriyoruz
    total_molar_mass = molar_mass_per_base * sequence_length  # DNA/RNA'nın toplam molar kütlesi
    copy_number_per_ul = (ng_in_grams * avogadro_number) / total_molar_mass  # Kopya sayısı hesaplama
    return copy_number_per_ul

# Seyreltme oranını hesaplama fonksiyonu
def calculate_dilution_factor(initial_cp_ul, target_cp_ul):
    if target_cp_ul >= initial_cp_ul:
        return 1, "Seyreltme gerekmiyor."
    dilution_factor = initial_cp_ul / target_cp_ul
    return dilution_factor, f"Önerilen seyreltme oranı: **1:{dilution_factor:.1f}**"

# Seri dilüsyon önerme fonksiyonu
def suggest_serial_dilution(dilution_factor):
    dilutions = [2, 5, 10]
    suggested = None
    for d in dilutions:
        steps = np.log(dilution_factor) / np.log(d)
        if steps.is_integer() or (steps % 1 < 0.2 or steps % 1 > 0.8):  # Yaklaşık tam sayı kontrolü
            suggested = f"{int(round(steps))} adımda **{d}-kat** seri dilüsyon yapabilirsiniz."
            break
    return suggested or "Seri dilüsyon için kesin bir oran bulunamadı, hassas seyreltme önerilir."

# Streamlit uygulaması başlığı
st.title("🔬 DNA/RNA Kopya Sayısı ve Seyreltme Hesaplayıcı")

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

    # Bilimsel gösterimde ve tam sayı olarak ayrıştırılmış gösterim
    cp_sci_notation = f"{cp_per_ul:.3e}"  # Bilimsel gösterim
    cp_digits = f"{cp_per_ul:,.0f}"  # Virgüllü tam sayı gösterimi

    # Sonucu çerçeve içinde gösterme
    st.markdown(
        f"""
        <div style="border: 2px solid #4CAF50; padding: 10px; border-radius: 10px; background-color: #e8f5e9; text-align: center; font-size: 20px;">
            <b>{ng_per_ul:.2f} ng/µL {strand_type}{molecule_type} için kopya sayısı:</b> <br>
            <b style="color: #2E7D32;">{cp_sci_notation} kopya/µL</b> <br>
            <b style="color: #1E88E5;">({cp_digits} kopya/µL)</b>
        </div>
        """, unsafe_allow_html=True
    )

    # Kullanıcıdan çalışmak istediği hedef kopya sayısını al (e formatı olmadan)
    target_cp_ul = st.number_input("Çalışmak istediğiniz kopya sayısını girin (cp/µL)", min_value=0.0, format="%.0f")

    if target_cp_ul > 0:
        dilution_factor, dilution_message = calculate_dilution_factor(cp_per_ul, target_cp_ul)
        serial_dilution_suggestion = suggest_serial_dilution(dilution_factor)

        # Seyreltme önerisini çerçeve içinde göster
        st.markdown(
            f"""
            <div style="border: 2px solid #FF9800; padding: 10px; border-radius: 10px; background-color: #FFF3E0; text-align: center; font-size: 18px;">
                <b>Hedef Kopya Sayısı: {target_cp_ul:,.0f} cp/µL</b> <br>
                <b style="color: #E65100;">{dilution_message}</b>
            </div>
            """, unsafe_allow_html=True
        )

        # Seri dilüsyon önerisini göster
        st.write(f"🔬 {serial_dilution_suggestion}")

else:
    st.write("Lütfen geçerli bir baz uzunluğu ve ng/µL değeri girin.")

# Hesaplama formülünü bilimsel olarak gösterme
st.subheader("📌 Hesaplama Formülü:")
st.latex(r"""
\text{Kopya Sayısı (cp/µL)} = \frac{\left( \text{Ng/µL} \times 10^{-9} \right) \times 6.022 \times 10^{23}}{\text{Baz Uzunluğu} \times \text{Molar Kütle (g/mol)}}
""")

st.write("""
### Açıklamalar:
- **Ng/µL**: Numune konsantrasyonu (nanogram/µL)
- **Avogadro Sayısı**: \\(6.022 \times 10^{23}\\) (Bir moldaki molekül sayısı)
- **Molar Kütle**: DNA/RNA’nın **bir baz başına** kütlesi (g/mol)  
  - ssDNA: **330 g/mol**  
  - dsDNA: **660 g/mol**  
  - ssRNA: **340 g/mol**  
- **Baz Uzunluğu**: DNA/RNA dizisinin uzunluğu (baz sayısı)
- **Seyreltme Oranı**: Başlangıç kopya sayısını hedeflenen kopya sayısına düşürmek için gerekli oran
- **Seri Dilüsyon**: Küçük adımlarla doğruluk sağlamak için önerilen ardışık seyreltme oranı
""")
