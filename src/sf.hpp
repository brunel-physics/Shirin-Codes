#ifndef SF_HPP
#define SF_HPP

#include <algorithm>
#include <array>

namespace sf
{
[[gnu::const]] inline float muon_id(const float pt, const float eta)
{
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017
    // https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ID.json
    // NUM_TightID_DEN_genTracks
    static constexpr std::array<float, 3> etabinedges{0.9f, 1.2f, 2.1f};
    static constexpr std::array<float, 5> ptbinedges{25.f, 30.f, 40.f, 50.f, 60.f};
    static constexpr std::array<std::array<std::pair<float, float>, 6>, 4> muon_id_sfs{{
        {{{0.9910777627756951f, 0.0034967203087024274f}, {0.9874104682620840f, 0.0018975082960536634f}, {0.9907753279135898f, 0.0003977503197109705f}, {0.9892483588952047f, 0.00032329941312374114f}, {0.9855545160334763f, 0.0008602730194340781f}, {0.9898057377093389f, 0.0016204327419630327f}}}, // η 0-0.9
        {{{0.9927389275515244f, 0.0047437370661283660f}, {0.9850639397625120f, 0.0218787608424192500f}, {0.9865359464182247f, 0.0006903526352667042f}, {0.9849130931014930f, 0.02013503091494561000f}, {0.9839056384760008f, 0.0015917233692683600f}, {0.9840604031434680f, 0.0121278780491192190f}}}, // η 0.9-1.2
        {{{0.9924252719877384f, 0.0077859527402420020f}, {0.9890884461284933f, 0.0149053560477533670f}, {0.9946469069883841f, 0.0124242236899199730f}, {0.9926528825155183f, 0.00993167811549696700f}, {0.9906364222943529f, 0.0009713213798502730f}, {0.9920464322143979f, 0.0021353964567237746f}}}, // η 1.2-2.1111
        {{{0.9758095839531763f, 0.0043993151217841040f}, {0.9745153594179884f, 0.0027111009825340473f}, {0.9787410500158746f, 0.0010035577872014160f}, {0.9781891229195010f, 0.00112306059413853970f}, {0.9673568416097894f, 0.0037006525169638958f}, {0.9766311856731202f, 0.0086266646688285500f}}}, // η 2.1+
    }};

    const auto eta_bin{std::distance(
            etabinedges.begin(),
            std::upper_bound(etabinedges.begin(), etabinedges.end(), eta))};
    const auto pt_bin{std::distance(
            ptbinedges.begin(),
            std::upper_bound(ptbinedges.begin(), ptbinedges.end(), pt))};

    return muon_id_sfs.at(eta_bin).at(pt_bin).first;
}

[[gnu::const]] inline float muon_iso(const float pt, const float eta)
{
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017
    // https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ISO.json
    // NUM_LooseRelIso_DEN_TightIDandIPCut
    static constexpr std::array<float, 3> etabinedges{0.9f, 1.2f, 2.1f};
    static constexpr std::array<float, 5> ptbinedges{25.f, 30.f, 40.f, 50.f, 60.f};
    static constexpr std::array<std::array<std::pair<float, float>, 6>, 4> muon_iso_sfs{{
        {{{0.9931118153898082f, 0.0023961047801277897f}, {0.9963959301142516f, 0.0011130262628559282f}, {0.9983181988469478f, 0.00028317863248081760f}, {0.9994347636372417f, 0.00011232714693801277f}, {0.9997680258722230f, 0.00021161394350289537f}, {1.0002119715073154f, 0.00034548765287024270f}}}, // η 0-0.9
        {{{0.9956674008369976f, 0.0040757770757467080f}, {0.9932947173775393f, 0.0020402913397330540f}, {0.9976914943066302f, 0.00052908143248948090f}, {0.9990342041383322f, 0.00020372497501125323f}, {0.9994155529201356f, 0.00041662994217661376f}, {1.0002046210970306f, 0.00066202613921456840f}}}, // η 0.9-1.2
        {{{0.9967414270112102f, 0.0016604845648881112f}, {0.9988489100332861f, 0.0009332375024769512f}, {0.9988826842051872f, 0.00028173209377228390f}, {0.9993872934568914f, 0.00011049799288742951f}, {0.9997410091519127f, 0.00024619498190323770f}, {1.0004725402685148f, 0.00041307624612303565f}}}, // η 1.2-2.1
        {{{0.9973954140213298f, 0.0021903221583127706f}, {0.9987560726170641f, 0.0012087380293472640f}, {0.9990618844158636f, 0.00037710084052130214f}, {0.9997635755214144f, 0.00017351798487330648f}, {1.0002224795137067f, 0.00051831556143466900f}, {0.9993431865975091f, 0.00085142066151194850f}}}, // η 2.1+
    }};

    const auto eta_bin{std::distance(
            etabinedges.begin(),
            std::upper_bound(etabinedges.begin(), etabinedges.end(), eta))};
    const auto pt_bin{std::distance(
            ptbinedges.begin(),
            std::upper_bound(ptbinedges.begin(), ptbinedges.end(), pt))};

    return muon_iso_sfs.at(eta_bin).at(pt_bin).first;
}
}

#endif /* SF_HPP */
