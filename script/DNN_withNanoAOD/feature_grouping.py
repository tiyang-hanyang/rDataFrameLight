from signalExtractionInputFromNanoAOD import FEATURE_FIELD_NAMES


GROUP_NAME_ALIASES = {
    "1": "muon_four_momentum",
    "2": "jet_four_momentum",
    "3": "jet_pair_dr_invmass",
    "4": "dimuon_dr_mass",
    "5": "muon_jet_dr",
    "6": "extra_muon_info",
    "7": "extra_jet_info",
}


def build_feature_groups():
    feature_indices = {name: index for index, name in enumerate(FEATURE_FIELD_NAMES)}
    groups = [
        (
            "muon_four_momentum",
            [
                feature_indices["lead_mu_px"],
                feature_indices["lead_mu_py"],
                feature_indices["lead_mu_pz"],
                feature_indices["sublead_mu_px"],
                feature_indices["sublead_mu_py"],
                feature_indices["sublead_mu_pz"],
            ],
        ),
        (
            "jet_four_momentum",
            [
                feature_indices[f"jet{jet_index}_{component}"]
                for jet_index in range(1, 5)
                for component in ("px", "py", "pz", "mass")
            ],
        ),
        (
            "jet_pair_dr_invmass",
            [
                feature_indices[f"bjet_dr_{pair[0]}{pair[1]}"]
                for pair in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
            ]
            + [
                feature_indices[f"bjet_m_{pair[0]}{pair[1]}"]
                for pair in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
            ],
        ),
        (
            "dimuon_dr_mass",
            [
                feature_indices["dimuon_dr"],
                feature_indices["dimuon_mass"],
            ],
        ),
        (
            "muon_jet_dr",
            [
                feature_indices[f"jet{jet_index}_dr_mu{mu_index}"]
                for jet_index in range(1, 5)
                for mu_index in (1, 2)
            ],
        ),
        (
            "extra_muon_info",
            [
                feature_indices["lead_mu_miniPFRelIso_all"],
                feature_indices["sublead_mu_miniPFRelIso_all"],
                feature_indices["lead_mu_jetRelIso"],
                feature_indices["sublead_mu_jetRelIso"],
                feature_indices["lead_mu_jetDF"],
                feature_indices["sublead_mu_jetDF"],
                feature_indices["lead_mu_promptMVA"],
                feature_indices["sublead_mu_promptMVA"],
            ],
        ),
        (
            "extra_jet_info",
            [
                feature_indices[f"jet{jet_index}_btag"]
                for jet_index in range(1, 5)
            ]
            + [
                feature_indices["nGoodJet"],
                feature_indices["nBJet"],
                feature_indices["sumjet_pt"],
                feature_indices["sumbjet_pt"],
                feature_indices["jet_cent"],
            ],
        ),
    ]

    assigned = sorted(index for _, indices in groups for index in indices)
    expected = list(range(len(FEATURE_FIELD_NAMES)))
    if assigned != expected:
        missing = [FEATURE_FIELD_NAMES[index] for index in expected if index not in assigned]
        duplicated = []
        seen = set()
        for _, indices in groups:
            for index in indices:
                if index in seen and FEATURE_FIELD_NAMES[index] not in duplicated:
                    duplicated.append(FEATURE_FIELD_NAMES[index])
                seen.add(index)
        raise RuntimeError(
            "Feature groups do not partition all inputs exactly once. "
            f"Missing={missing}, duplicated={duplicated}"
        )
    return groups


def get_group_names():
    return [group_name for group_name, _ in build_feature_groups()]


def normalize_group_token(token):
    token = str(token).strip()
    if token in GROUP_NAME_ALIASES:
        return GROUP_NAME_ALIASES[token]
    return token


def resolve_feature_groups(group_tokens):
    all_group_names = get_group_names()
    if group_tokens is None:
        return list(all_group_names)

    resolved = []
    for token in group_tokens:
        group_name = normalize_group_token(token)
        if group_name not in all_group_names:
            raise RuntimeError(
                f"Unknown feature group '{token}'. Available: {', '.join(all_group_names)} "
                "or numeric aliases 1..7."
            )
        if group_name not in resolved:
            resolved.append(group_name)

    if not resolved:
        raise RuntimeError("At least one feature group must be selected.")
    return resolved


def build_feature_mask(active_group_names):
    groups = build_feature_groups()
    active = set(resolve_feature_groups(active_group_names))
    mask = [False] * len(FEATURE_FIELD_NAMES)
    for group_name, feature_indices in groups:
        if group_name not in active:
            continue
        for feature_index in feature_indices:
            mask[feature_index] = True
    return mask


def get_active_feature_indices(active_group_names):
    groups = build_feature_groups()
    active = set(resolve_feature_groups(active_group_names))
    indices = []
    for group_name, feature_indices in groups:
        if group_name in active:
            indices.extend(feature_indices)
    return sorted(indices)


def get_active_feature_names(active_group_names):
    return [FEATURE_FIELD_NAMES[index] for index in get_active_feature_indices(active_group_names)]


def apply_feature_mask(features, active_group_names):
    resolved_groups = resolve_feature_groups(active_group_names)
    mask = build_feature_mask(resolved_groups)
    masked = features.astype(features.dtype, copy=True)
    for feature_index, keep in enumerate(mask):
        if not keep:
            masked[:, feature_index] = 0.0
    return masked


def select_active_features(features, active_group_names):
    indices = get_active_feature_indices(active_group_names)
    return features[:, indices].astype(features.dtype, copy=False)


def describe_feature_groups(active_group_names):
    active = set(resolve_feature_groups(active_group_names))
    descriptions = []
    for group_name, feature_indices in build_feature_groups():
        if group_name in active:
            descriptions.append(
                {
                    "group_name": group_name,
                    "n_features": len(feature_indices),
                    "features": [FEATURE_FIELD_NAMES[index] for index in feature_indices],
                }
            )
    return descriptions
