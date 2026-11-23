# b por macroregión, tomado de la Extended Data Table 4 del paper
# (Bastien-Olvera et al., Nature 2024)
using CSV, DataFrames

const REGION_B = Dict(
    "OECD" => 0.153859,
    "LAM"  => 0.193776,
    "MAF"  => 0.322886,
    "ASIA" => 0.559234,
    "REF"  => 0.612428,
)

const REGION_COUNTRIES = Dict(
    "OECD" => ["ALB","AUS","AUT","BEL","BGR","BIH","CAN","CHE","CHL","COL","CRI",
               "CYP","CZE","DEU","DNK","ESP","EST","FIN","FRA","GBR","GRC","HRV",
               "HUN","IRL","ISL","ISR","ITA","JPN","KOR","LTU","LUX","LVA","MDA",
               "MEX","MKD","MLT","MNE","NLD","NOR","NZL","POL","PRT","ROU","SRB",
               "SVK","SVN","SWE","TUR","UKR","USA"],
    "LAM"  => ["ABW","ARG","BHS","BLZ","BOL","BRA","BRB","CUB","DOM","ECU","GTM",
               "GUY","HND","HTI","JAM","LCA","NIC","PAN","PER","PRI","PRY","SLV",
               "SUR","TTO","URY","VCT","VEN"],
    "ASIA" => ["AFG","BGD","BRN","BTN","CHN","FJI","HKG","IDN","IND","KHM","LAO",
               "LKA","MAC","MDV","MMR","MNG","MYS","NCL","NPL","PAK","PHL","PNG",
               "PYF","SGP","SLB","THA","TLS","TON","TWN","VNM","VUT","WSM"],
    "REF"  => ["ARM","AZE","BLR","GEO","KAZ","KGZ","RUS","TJK","TKM","UZB"],
    "MAF"  => ["AGO","ARE","BDI","BEN","BFA","BHR","BWA","CAF","CIV","CMR","COD",
               "COG","COM","CPV","DJI","DZA","EGY","ERI","ETH","GAB","GHA","GIN",
               "GMB","GNB","GNQ","IRN","IRQ","JOR","KEN","KWT","LBN","LBR","LBY",
               "LSO","MAR","MDG","MLI","MOZ","MRT","MUS","MWI","NAM","NER","NGA",
               "OMN","PSE","QAT","RWA","SAU","SDN","SEN","SLE","SOM","STP","SWZ",
               "SYR","TCD","TGO","TUN","TZA","UGA","YEM","ZAF","ZMB","ZWE"],
)

# País → macroregión, con detección de duplicados
const COUNTRY2REGION = Dict{String,String}()
for (r, isos) in REGION_COUNTRIES, iso in isos
    haskey(COUNTRY2REGION, iso) && error("$iso asignado a dos regiones")
    COUNTRY2REGION[iso] = r
end

# Vector de b en el mismo orden que `countries` de get_model
function b_vector(root::AbstractString)
    countries = String.(CSV.read(joinpath(root, "data", "Dimension_countries.csv"),
                                 DataFrame).CountryISO)
    missing_iso = setdiff(countries, keys(COUNTRY2REGION))
    extra_iso   = setdiff(keys(COUNTRY2REGION), countries)
    isempty(missing_iso) || error("Sin región: $missing_iso")
    isempty(extra_iso)   || @warn "Códigos que no están en el modelo: $extra_iso"
    return [REGION_B[COUNTRY2REGION[c]] for c in countries]
end