Circos <- S7::new_class(
  "Circos",
  properties = list(

    # Writable
    zones = S7::class_list,
    inter_chrom_spacing = S7::class_numeric,

    # Computed
    n_zones = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) { length(self@zones)},
      setter = function(self, value) { "@n_zones is a read-only property" }
    ),
    zonetypes = S7::new_property(
      class = S7::class_character,
      getter = function(self) {
        vapply(self@zones, FUN = function(z) { z@type}, FUN.VALUE = character(1))
      },
      setter = function(self, value) { "@zonetypes is a read-only property" }
    ),
    n_reference_zones = S7::new_property(
      class = S7::class_numeric,
      getter = function(self) { sum(self@zonetypes == "Reference") },
      setter = function(self, value) { "@n_reference_zones is a read-only property" }
    ),

    # Functions
    project_point = S7::new_property(
      class = S7::class_function,
      getter = function(self){
        self@n_reference_zones
      }
    )
  ),
  validator = function(self) {
    if (self@n_reference_zones != 1) sprintf("Circos class object must have exactly one reference genome, not %i", self@n_reference_zones)
    if (self@zonetypes[[1]] != "Reference") sprintf("First zone must be of type `Reference`, not %s", self@zonetypes[[1]])
    zones_classes_valid <- all(vapply(self@zones, function(z) {inherits(z, "ggcircos::Zone") }, FUN.VALUE = logical(1)))
    if (!zones_classes_valid) sprintf("@zones paramater must be entirely composed of ggcircos::Zone objects created by `Zone()` function")
  },
  constructor = function(reference, inter_chrom_spacing = 0) {

    # Identify reference
    reference_zone <- Zone(reference, type = "Reference")

    S7::new_object(
      S7::S7_object(),
      zones = list(reference_zone),
      inter_chrom_spacing = inter_chrom_spacing
    )
})

Zone <- S7::new_class("Zone",
  properties = list(
    type = S7::new_property(
      class = S7::class_character,
      validator = function(value) {
        if (length(value) != 1) {
          return("@type must be a string (length 1 character vector)")
        }
        if (!value %in% valid_zone_types()) {
          return(sprintf("@type must be one of [%s]", toString(valid_zone_types())))
        }
      }
    ),
    # A data.frame with columns whose expected cols depend on type
    data = S7::new_property(
      class = S7::class_data.frame
    )
  ),
  constructor = function(data, type) {
    assertions::assert_one_of(type, valid_zone_types())
    assertions::assert_dataframe(data)

    if (type == "Reference") {
      assertions::assert_names_include(data, c("name", "start", "end", "refindex"))
    }

    if (type == "SmallMutations") {
      assertions::assert_names_include(data, c("name", "position", "vaf"))
    }

    S7::new_object(
      S7::S7_object(),
      data = data,
      type = type
    )
  }
)

valid_zone_types <- function() {
  c("Reference", "Copynumber", "SmallMutations", "Links")
}
