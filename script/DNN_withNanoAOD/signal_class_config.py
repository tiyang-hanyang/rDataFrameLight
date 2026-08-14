from dataclasses import dataclass
from typing import Optional, Tuple


TTHH_DL_ALIASES = ("TTHH_DL", "TTHH_2B2W_DL")
TTHH_SL_ALIASES = ("TTHH_SL", "TTHH_2B2W_SL")


@dataclass(frozen=True)
class ChannelSpec:
    channel_name: str
    class_name: str
    aliases: Tuple[str, ...]
    balance_group: Optional[str] = None


@dataclass(frozen=True)
class SignalClassSchema:
    name: str
    description: str
    class_names: Tuple[str, ...]
    channels: Tuple[ChannelSpec, ...]
    sample_rules: Tuple[Tuple[str, str], ...]
    allowed_os_channels: Optional[Tuple[str, ...]] = None
    allowed_sr_channels: Optional[Tuple[str, ...]] = None

    @property
    def class_to_index(self):
        return {class_name: index for index, class_name in enumerate(self.class_names)}

    @property
    def channel_to_class(self):
        return {channel.channel_name: channel.class_name for channel in self.channels}

    def infer_label(self, sample_name):
        sample_upper = sample_name.upper()
        for token, class_name in self.sample_rules:
            if token in sample_upper:
                return self.class_to_index[class_name]
        raise RuntimeError(f"Unable to infer label from sample name: {sample_name}")

    def identify_channel(self, name):
        upper_name = str(name).upper()
        for channel in self.channels:
            if any(alias in upper_name for alias in channel.aliases):
                return channel.channel_name, channel.class_name
        raise RuntimeError(f"Cannot determine channel from name: {name}")

    def canonicalize_channel_name(self, name):
        channel_name, _ = self.identify_channel(name)
        return channel_name

    def is_allowed_in_os(self, channel_name):
        if self.allowed_os_channels is None:
            return True
        return channel_name in self.allowed_os_channels

    def is_allowed_in_sr(self, channel_name):
        if self.allowed_sr_channels is None:
            return True
        return channel_name in self.allowed_sr_channels


THREE_CLASS_SCHEMA = SignalClassSchema(
    name="three_class",
    description="3 classes: TTHH | ttbar+TTBB | TTW/TTZ/TTTT/TTH",
    class_names=("TTHH", "ttbar_ttbb", "ttX_like"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("ttbarDL", "ttbar_ttbb", ("TTBARDL", "TTBAR_DL")),
        ChannelSpec("ttbarSL", "ttbar_ttbb", ("TTBARSL", "TTBAR_SL")),
        ChannelSpec("TTBB_DL", "ttbar_ttbb", ("TTBB_DL",)),
        ChannelSpec("TTBB_SL", "ttbar_ttbb", ("TTBB_SL",)),
        ChannelSpec("TTW", "ttX_like", ("TTW",)),
        ChannelSpec("TTZ_low", "ttX_like", ("TTZ_LOW", "TTZLOW")),
        ChannelSpec("TTZ_high", "ttX_like", ("TTZ_HIGH", "TTZHIGH")),
        ChannelSpec("TTTT", "ttX_like", ("TTTT",)),
        ChannelSpec("TTHBB", "ttX_like", ("TTHBB",)),
        ChannelSpec("TTHnonBB", "ttX_like", ("TTHNONBB", "TTH_NONBB")),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTBAR", "ttbar_ttbb"),
        ("TTBB", "ttbar_ttbb"),
        ("TTW", "ttX_like"),
        ("TTZ", "ttX_like"),
        ("TTTT", "ttX_like"),
        ("TTH", "ttX_like"),
        ("THQ", "ttX_like"),
        ("THW", "ttX_like"),
    ),
)

THREE_CLASS_OMIT_TTZ_MERGED_TTTT_SCHEMA = SignalClassSchema(
    name="three_class_omit_ttz_merged_tttt",
    description=(
        "3 classes: TTHH | tt_b(ttbar+TTBB+TTHBB) | ttX_like(TTW+TTHnonBB+TTTT), "
        "with TTZ removed from all OS splits and kept only in SRtest"
    ),
    class_names=("TTHH", "tt_b", "ttX_like"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("ttbar", "tt_b", ("TTBARDL", "TTBAR_DL", "TTBARSL", "TTBAR_SL")),
        ChannelSpec("TTBB", "tt_b", ("TTBB_DL", "TTBB_SL")),
        ChannelSpec("TTHBB", "tt_b", ("TTHBB",)),
        ChannelSpec("TTW", "ttX_like", ("TTW",), "ttX_like_merged"),
        ChannelSpec("TTHnonBB", "ttX_like", ("TTHNONBB", "TTH_NONBB"), "ttX_like_merged"),
        ChannelSpec("TTTT", "ttX_like", ("TTTT",), "ttX_like_merged"),
        ChannelSpec("TTZ", "ttX_like", ("TTZ_LOW", "TTZLOW", "TTZ_HIGH", "TTZHIGH", "TTZ"), "TTZ"),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTHNONBB", "ttX_like"),
        ("TTHBB", "tt_b"),
        ("TTBAR", "tt_b"),
        ("TTBB", "tt_b"),
        ("TTTT", "ttX_like"),
        ("TTW", "ttX_like"),
        ("TTZ", "ttX_like"),
        ("TTH", "ttX_like"),
        ("THQ", "ttX_like"),
        ("THW", "ttX_like"),
    ),
    allowed_os_channels=("TTHH_DL", "TTHH_SL", "ttbar", "TTBB", "TTHBB", "TTW", "TTHnonBB", "TTTT"),
    allowed_sr_channels=("TTHH_DL", "TTHH_SL", "ttbar", "TTBB", "TTHBB", "TTW", "TTZ", "TTHnonBB", "TTTT"),
)

BINARY_TTHH_VS_TTB_SCHEMA = SignalClassSchema(
    name="binary_tthh_vs_ttb",
    description="2 classes: TTHH | ttbar+TTBB+TTHBB",
    class_names=("TTHH", "tt_b"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("ttbar", "tt_b", ("TTBARDL", "TTBAR_DL", "TTBARSL", "TTBAR_SL")),
        ChannelSpec("TTBB", "tt_b", ("TTBB_DL", "TTBB_SL")),
        ChannelSpec("TTHBB", "tt_b", ("TTHBB",)),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTHBB", "tt_b"),
        ("TTBAR", "tt_b"),
        ("TTBB", "tt_b"),
    ),
)

BINARY_TTHH_VS_TTWLIKE_SCHEMA = SignalClassSchema(
    name="binary_tthh_vs_ttw_tthnonbb",
    description="2 classes: TTHH | TTW+TTHnonBB",
    class_names=("TTHH", "ttW_like"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("TTW", "ttW_like", ("TTW",)),
        ChannelSpec("TTHnonBB", "ttW_like", ("TTHNONBB", "TTH_NONBB")),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTHNONBB", "ttW_like"),
        ("TTW", "ttW_like"),
    ),
)

BINARY_TTHH_VS_TTTT_SCHEMA = SignalClassSchema(
    name="binary_tthh_vs_tttt",
    description="2 classes: TTHH | TTTT",
    class_names=("TTHH", "TTTT"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("TTTT", "TTTT", ("TTTT",)),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTTT", "TTTT"),
    ),
)

FOUR_CLASS_MERGED_SCHEMA = SignalClassSchema(
    name="four_class_merged",
    description="4 classes: TTHH | ttbar+TTBB+TTHBB | TTW+TTZ+TTHnonBB | TTTT, with TTZ_low/high merged into TTZ",
    class_names=("TTHH", "tt_b", "ttX_like", "TTTT"),
    channels=(
        ChannelSpec("TTHH_DL", "TTHH", TTHH_DL_ALIASES),
        ChannelSpec("TTHH_SL", "TTHH", TTHH_SL_ALIASES),
        ChannelSpec("ttbar", "tt_b", ("TTBARDL", "TTBAR_DL", "TTBARSL", "TTBAR_SL")),
        ChannelSpec("TTBB", "tt_b", ("TTBB_DL", "TTBB_SL")),
        ChannelSpec("TTHBB", "tt_b", ("TTHBB",)),
        ChannelSpec("TTW", "ttX_like", ("TTW",)),
        ChannelSpec("TTZ", "ttX_like", ("TTZ_LOW", "TTZLOW", "TTZ_HIGH", "TTZHIGH", "TTZ")),
        ChannelSpec("TTHnonBB", "ttX_like", ("TTHNONBB", "TTH_NONBB")),
        ChannelSpec("TTTT", "TTTT", ("TTTT",)),
    ),
    sample_rules=(
        ("TTHH", "TTHH"),
        ("TTHNONBB", "ttX_like"),
        ("TTHBB", "tt_b"),
        ("TTBAR", "tt_b"),
        ("TTBB", "tt_b"),
        ("TTTT", "TTTT"),
        ("TTW", "ttX_like"),
        ("TTZ", "ttX_like"),
        ("TTH", "ttX_like"),
        ("THQ", "ttX_like"),
        ("THW", "ttX_like"),
    ),
)

FOUR_CLASS_OMIT_TTZ_SCHEMA = SignalClassSchema(
    name="four_class_omit_ttz",
    description=(
        "4 classes: TTHH | ttbar+TTBB+TTHBB | TTW+TTZ+TTHnonBB | TTTT, "
        "but TTZ is removed from all OS splits and kept only in SRtest"
    ),
    class_names=("TTHH", "tt_b", "ttX_like", "TTTT"),
    channels=FOUR_CLASS_MERGED_SCHEMA.channels,
    sample_rules=FOUR_CLASS_MERGED_SCHEMA.sample_rules,
    allowed_os_channels=("TTHH_DL", "TTHH_SL", "ttbar", "TTBB", "TTHBB", "TTW", "TTHnonBB", "TTTT"),
    allowed_sr_channels=("TTHH_DL", "TTHH_SL", "ttbar", "TTBB", "TTHBB", "TTW", "TTZ", "TTHnonBB", "TTTT"),
)


SCHEMAS = {
    THREE_CLASS_SCHEMA.name: THREE_CLASS_SCHEMA,
    THREE_CLASS_OMIT_TTZ_MERGED_TTTT_SCHEMA.name: THREE_CLASS_OMIT_TTZ_MERGED_TTTT_SCHEMA,
    BINARY_TTHH_VS_TTB_SCHEMA.name: BINARY_TTHH_VS_TTB_SCHEMA,
    BINARY_TTHH_VS_TTWLIKE_SCHEMA.name: BINARY_TTHH_VS_TTWLIKE_SCHEMA,
    BINARY_TTHH_VS_TTTT_SCHEMA.name: BINARY_TTHH_VS_TTTT_SCHEMA,
    FOUR_CLASS_MERGED_SCHEMA.name: FOUR_CLASS_MERGED_SCHEMA,
    FOUR_CLASS_OMIT_TTZ_SCHEMA.name: FOUR_CLASS_OMIT_TTZ_SCHEMA,
}

DEFAULT_SCHEMA_NAME = THREE_CLASS_OMIT_TTZ_MERGED_TTTT_SCHEMA.name

CATEGORY_COLORS = {
    "TTHH": "#d62728",
    "ttbar_ttbb": "#1f77b4",
    "ttbar_like": "#1f77b4",
    "tt_b": "#1f77b4",
    "ttW_like": "#2ca02c",
    "ttX_like": "#2ca02c",
    "TTTT": "#ff7f0e",
}

CHANNEL_COLORS = {
    "TTHH_DL": "#d62728",
    "TTHH_SL": "#ff9896",
    "ttbarDL": "#1f77b4",
    "ttbarSL": "#9ecae1",
    "ttbar": "#1f77b4",
    "TTBB_DL": "#17becf",
    "TTBB_SL": "#9edae5",
    "TTBB": "#17becf",
    "TTW": "#2ca02c",
    "TTZ_low": "#98df8a",
    "TTZ_high": "#006d2c",
    "TTZ": "#2ca25f",
    "TTTT": "#ff7f0e",
    "TTHBB": "#9467bd",
    "TTHnonBB": "#8c564b",
}


def get_signal_schema(name=None):
    resolved_name = DEFAULT_SCHEMA_NAME if name is None else name
    if resolved_name not in SCHEMAS:
        raise RuntimeError(
            f"Unknown class schema: {resolved_name}. Available: {', '.join(sorted(SCHEMAS))}"
        )
    return SCHEMAS[resolved_name]


def available_schema_names():
    return tuple(sorted(SCHEMAS))


def schema_help_text():
    parts = [f"{schema.name} ({schema.description})" for schema in SCHEMAS.values()]
    return "Signal class schema. Available: " + "; ".join(parts)


def get_category_color(name):
    return CATEGORY_COLORS.get(name, "#4c4c4c")


def get_channel_color(name):
    return CHANNEL_COLORS.get(name, "#4c4c4c")
