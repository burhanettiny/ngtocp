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
        return 1, "No dilution required."
    dilution_factor = initial_cp_ul / target_cp_ul
    return dilution_factor, f"Suggested dilution factor: **1:{dilution_factor:.1f}**"

# Digital PCR reaksiyonu için kopya sayısı hesaplama
def calculate_digital_pcr_copies(cp_per_ul, template_volume_ul, reaction_volume_ul):
    copies_per_reaction = cp_per_ul * template_volume_ul  # Reaksiyona giren toplam kopya
    copies_per_ul_reaction = copies_per_reaction / reaction_volume_ul  # Reaksiyon hacmine normalize etme
    return copies_per_reaction, copies_per_ul_reaction

# Dil seçeneği
language = st.radio("Dil Seçin:", ("Türkçe", "English"))

if language == "Türkçe":
    # Başlık ve metinler Türkçe
    st.title("🔬 DNA/RNA Kopya Sayısı Hesaplayıcı")
    molecule_type_label = "Molekül tipi:"
    sequence_length_label = "Baz uzunluğunu girin (baz sayısı)"
    ng_per_ul_label = "Konsantrasyonu girin (ng/µL)"
    target_cp_ul_label = "Çalışmak istediğiniz kopya sayısını girin (cp/µL)"
    reaction_volume_label = "Digital PCR reaksiyon hacmi (µL)"
    template_volume_label = "Template DNA/RNA hacmi (µL)"
    dilution_message_default = "Seyreltme gerekmiyor."
    dilution_message_suggestion = "Önerilen seyreltme oranı: **1:{:.1f}**"
    digital_pcr_message = "Digital PCR için toplam kopya sayısı:"
    formula_title = "📌 Hesaplama Formülü:"
    formula_1 = r"\text{Kopya Sayısı (cp/µL)} = \frac{\left( \text{Ng/µL} \times 10^{-9} \right) \times 6.022 \times 10^{23}}{\text{Baz Uzunluğu} \times \text{Molar Kütle (g/mol)}}"
    formula_2 = r"\text{Digital PCR Kopya Sayısı} = \text{(Kopya/µL)} \times \text{Template Hacmi (µL)}"
else:
    # Başlık ve metinler İngilizce
    st.title("🔬 DNA/RNA Copy Number Calculator")
    molecule_type_label = "Molecule type:"
    sequence_length_label = "Enter sequence length (base count)"
    ng_per_ul_label = "Enter concentration (ng/µL)"
    target_cp_ul_label = "Enter target copy number (cp/µL)"
    reaction_volume_label = "Digital PCR reaction volume (µL)"
    template_volume_label = "Template DNA/RNA volume (µL)"
    dilution_message_default = "No dilution required."
    dilution_message_suggestion = "Suggested dilution factor: **1:{:.1f}**"
    digital_pcr_message = "Total copy number for Digital PCR:"
    formula_title = "📌 Calculation Formula:"
    formula_1 = r"\text{Copy Number (cp/µL)} = \frac{\left( \text{Ng/µL} \times 10^{-9} \right) \times 6.022 \times 10^{23}}{\text{Sequence Length} \times \text{Molar Mass (g/mol)}}"
    formula_2 = r"\text{Digital PCR Copy Number} = \text{(Copy/µL)} \times \text{Template Volume (µL)}"

# Kullanıcıdan DNA veya RNA türünü seçmesini iste
molecule_type = st.radio(molecule_type_label, ("DNA", "RNA"))

# Eğer DNA seçildiyse ss veya ds seçeneğini göster
if molecule_type == "DNA":
    strand_type = st.radio("Is it single-stranded (ss) or double-stranded (ds)?", ("ss", "ds"))
    molar_mass_per_base = 660 if strand_type == "ds" else 330  # dsDNA = 660 g/mol/baz, ssDNA = 330 g/mol/baz
else:
    strand_type = "ss"  # RNA genellikle tek zincirli olduğundan otomatik olarak "ss" olarak atanır.
    molar_mass_per_base = 340  # ssRNA = 340 g/mol/baz

# Kullanıcıdan baz uzunluğu girmesini iste
sequence_length = st.number_input(sequence_length_label, min_value=1, value=1000)

# Kullanıcıdan ng/uL cinsinden miktar girmesini iste
ng_per_ul = st.number_input(ng_per_ul_label, min_value=0.0, format="%.2f")

# Hesaplama yap
if sequence_length > 0 and ng_per_ul > 0:
    cp_per_ul = ng_to_cp_per_ul(ng_per_ul, sequence_length, molar_mass_per_base)

    # Bilimsel gösterimde ve tam sayı olarak ayrıştırılmış gösterim
    cp_sci_notation = f"{cp_per_ul:.3e}"  # Bilimsel gösterim
    cp_digits = f"{cp_per_ul:.0f}"  # Virgülsüz tam sayı gösterimi

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

    # Digital PCR hesaplama bölümü
    st.subheader("🧬 Digital PCR Reaksiyon Hesaplaması")

    # Kullanıcıdan dPCR reaksiyon hacmi ve template hacmi al
    reaction_volume_ul = st.number_input(reaction_volume_label, min_value=0.1, value=20.0, step=0.1)
    template_volume_ul = st.number_input(template_volume_label, min_value=0.1, value=2.0, step=0.1)

    if reaction_volume_ul > 0 and template_volume_ul > 0:
        copies_per_reaction, copies_per_ul_reaction = calculate_digital_pcr_copies(cp_per_ul, template_volume_ul, reaction_volume_ul)

        # Sonucu çerçeve içinde gösterme
        st.markdown(
            f"""
            <div style="border: 2px solid #3F51B5; padding: 10px; border-radius: 10px; background-color: #E3F2FD; text-align: center; font-size: 20px;">
                <b>{digital_pcr_message}</b> <br>
                <b style="color: #1A237E;">{copies_per_reaction:.0f} kopya/reaksiyon</b> <br>
                <b style="color: #304FFE;">({copies_per_ul_reaction:.2f} kopya/µL.reaksiyon hacmi)</b>
              </div>
            """, unsafe_allow_html=True
        )
    # Kullanıcıdan çalışmak istediği hedef kopya sayısını al (e formatı olmadan)
    target_cp_ul = st.number_input(target_cp_ul_label, min_value=0.0, format="%.0f")

    if target_cp_ul > 0:
        dilution_factor, dilution_message = calculate_dilution_factor(cp_per_ul, target_cp_ul)

        # Seyreltme önerisini çerçeve içinde göster
        st.markdown(
            f"""
            <div style="border: 2px solid #FF9800; padding: 10px; border-radius: 10px; background-color: #FFF3E0; text-align: center; font-size: 18px;">
                <b>Hedef Kopya Sayısı: {target_cp_ul:.0f} cp/µL</b> <br>
                <b style="color: #E65100;">{dilution_message}</b>
            </div>
            """, unsafe_allow_html=True
        )

# Hesaplama formülünü bilimsel olarak gösterme
st.subheader(formula_title)
st.latex(formula_1)
st.latex(formula_2)
