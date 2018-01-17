source "http://rubygems.org"
# Add dependencies required to use your gem here.
# Example:
#   gem "activesupport", ">= 2.3.5"

gem "bio", ">= 1.5.1"
gem "bio-samtools", ">= 2.6.2"
#gem "rake"

gem "systemu", ">=2.5.2"

group :development do
	gem "shoulda", ">= 2.10"
	gem 'test-unit'
	if RUBY_VERSION.start_with?("2.1") or RUBY_VERSION.start_with?("2.2") or RUBY_VERSION.start_with?("2.0")
		gem "jeweler", "= 2.0.1"
	else
		gem "juwelier" 
	end
end